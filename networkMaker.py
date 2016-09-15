#!/usr/env/python -c
import argparse
import sys
sys.path.append('/usr/local/share/py/djscripts') ;
import os
from os import listdir
from io import StringIO
import colorama
import numpy as np
import yaml
import re
from lib import interactors as I
from lib import interactors_extras as ie
from lib import rbase
from lib.markutils import b4us,afterus
from lib import config as cf

DB           = False
deprecation  = colorama.Fore.RED+" DEPRECATED 21 APRIL 2016 "+colorama.Fore.RESET
config       = { 'ifiles' : cf.ifilesPath,
                 'publicDatadir' : cf.publicDataPath,
                 'youtfilename' : '',
                 'outfilename'  : '',
                 'rescue_f'     : None }

def tokey(c, s) :
    if c['organism'] == 'human':
        if s in rbase.hmg['Symbol'] : 
            return s+'_'+rbase.hmg['Symbol'][s]['EID']
        elif re.search( '^[a-zA-Z0-9]+_\d+$', s ):
            return s
        else : 
            return s+'_00' ;
    elif c['organism'] == 'mouse':
        if s in rbase.mmg['Symbol'] : 
            return s+'_'+rbase.mmg['Symbol'][s]['EID']
        elif re.search( '^[a-zA-Z0-9]+_\d+$', s ):
            return s
        else : 
            return s+'_00' ;

    
def loadObjects( c ):

    if c['organism'] == 'human' :
        rbase.load('hmg')
    elif c['organism'] == 'mouse' :
        rbase.load('mmg')

        
def readYAMLfile( yamlfile, c ) :

    with open( yamlfile ) as infile :
        yout          = yaml.load(infile.read())

    c['ds_dicts']     = yout['datasets']

    refFiles          = dict()
    for dsd in c['ds_dicts'] : 
        lost_files    = 0
        if not os.path.isfile('./'+dsd['infilename']) and \
           not os.path.isfile(c['ifiles']+dsd['infilename']) : 
            print('File '+dsd['infilename']+' not found.')
            lost_files += 1
        if dsd['control'] not in refFiles:
            pFile     = ie.createReferenceFile( dsd['control'], overWrite = False)
            refFiles[dsd['control']] = pFile
            dsd['control']           = pFile
                                                 
    if lost_files > 0 : 
        raise IOError

    c['organism']     = yout['options'].get('organism','human')
    c['ALPHA_HI']     = yout['options'].get('alpha_hi',0.01)
    c['ALPHA_LO']     = yout['options'].get('alpha_lo',0.05)
    c['FLOOR_HI']     = yout['options'].get('floor_hi',4)
    c['FLOOR_LO']     = yout['options'].get('floor_lo',4)
    c['CORRL_HI']     = yout['options'].get('corrl_hi',70)
    c['CORRL_LO']     = yout['options'].get('corrl_lo',70)
    c['rescue_deg']   = yout['options'].get('degree_to_rescue',1)
    c['valid_quals']  = yout['options'].get('valid_degree_quals',{'wt',})
    c['mt_method']    = yout['options'].get('mt_method','fdr_bh')
    c['node_filter']  = yout['options'].get('node_filter', None)
    c['iact_filter']  = yout['options'].get('iact_filter', None) 
    c['nwd']          = yout['options'].get('network-wide_degree',False)
    c['public_dicts'] = yout['public']
    # this should be a list of dicts
    # each dict (1/public file) should have the fields : infilename qualify convert misncore minweight
    # public datasets have NO baits, are NOT DIRECTED, and are NEVER compared to negative controls

    if type(c['valid_quals']) is list : 
        c['valid_quals']  = set(c['valid_quals'])

    if yout.get('files') : 
        c['rescue_f']     = yout['files'].get('rescue',None)
        c['youtfilename'] = yout['files'].get('outfile',None)
        c['yidbfilename'] = yout['files'].get('idb',None)

    if c['rescue_f'] :
        if not os.path.isfile( c['rescue_f'] ) : 
            sys.stderr.write('Invalid path provided for rescue file.\n')
            os._exit(1)
    if not 'youtfilename' in c :
        sys.stderr.write('No valid output file name provided!\n')
        os._exit(1)
    else:
        c['outfilename']  = c['youtfilename']

    if 'yidbfilename' in c :
        c['idbfilename']  = c['yidbfilename']
        

def readInDatasets( nwdata, c ):

    # parse datasets
    baitkeys      = list()

    # iterating through the dataset items
    for dsd in c['ds_dicts'] : 

        if dsd['infilename'] in os.listdir('.') : 
            dsf = open(dsd['infilename'])
        elif os.path.isfile( c['ifiles'] + dsd['infilename'] ): 
            dsf = open(c['ifiles'] + dsd['infilename'])
        else :
            # problem
            continue
        
        if dsd.get('convert') == 'm2h' : 
            nwdata.parse( dsf, fd = I.fdms, m2h = True, qualify = dsd.get('qualify',''), directed = True )
        elif dsd.get('convert') == 'h2m' : 
            nwdata.parse( dsf, fd = I.fdms, h2m = True, qualify = dsd.get('qualify',''), directed = True )
        else : 
            nwdata.parse( dsf, fd = I.fdms, qualify = dsd.get('qualify',''), directed = True )

        baitkeys.append(tokey(c, dsd['bait']))

        dsf.close()

    c['baitkeys'] = baitkeys
    

def filterNodesByBackground( nwdata, c ):

    zfhits_strong          = set( )
    zfhits_weak            = set( )
    
    for dsd in c['ds_dicts'] : 
        zfhits_strong     |= ie.madfilter_corr( nwdata, dsd['control'], tokey(c, dsd['bait']),
                                                qual = dsd.get('qualify'), directed = True, alpha = c['ALPHA_HI'],
                                                floor = c['FLOOR_HI'], maxcorr = c['CORRL_HI'], debug = DB , 
                                                method = c['mt_method'] )
        zfhits_weak       |= ie.madfilter_corr( nwdata, dsd['control'], tokey(c, dsd['bait']),
                                                qual = dsd.get('qualify'), directed = True, alpha = c['ALPHA_LO'],
                                                floor = c['FLOOR_LO'], maxcorr = c['CORRL_LO'], debug = DB ,
                                                method = c['mt_method'] )
        
    # expt bait to node edges, that are deemed significant
    zfhits_joint           = zfhits_strong | zfhits_weak

    # pass 1 nodes : every node at the end of one of the validated edges
    # remove nodes that don't have experimental edges pointing to them
    node_pass1_all         = { nk for ek in zfhits_joint  for nk in { nwdata.edges[ek].to.key, nwdata.edges[ek].whence.key}}
    node_pass1_strong      = { nk for ek in zfhits_strong for nk in { nwdata.edges[ek].to.key, nwdata.edges[ek].whence.key}}

    c['zfhits_joint']      = zfhits_joint
    c['node_pass1_all']    = node_pass1_all
    c['node_pass1_strong'] = node_pass1_strong
    

def readPublicDatasets( nwdata, c ):

    # OK, NOW we dump in the public datasets
    for pd in c['public_dicts'] :
        # (done) make emili follow biogrid field conventions
        # TODO ditto bioplex
        if pd['infilename'] in os.listdir('.') : 
            pdsf = open(pd['infilename'])
        elif os.path.isfile(pd['infilename']) :
            pdsf = open(pd['infilename'])
        else : 
            pdsf = open(c['publicDatadir'] + pd['infilename'])

        temporaryds = I.dataSet( i_filter = c['iact_filter'], debug = DB )
        if pd.get('convert') == 'm2h' : 
            temporaryds.parse( pdsf, fd = I.fd_biogrid, m2h = True, qualify = pd.get('qualify',''),
                               directed = False, force_qualify = True )
        elif pd.get('convert') == 'h2m' : 
            temporaryds.parse( pdsf, fd = I.fd_biogrid, h2m = True, qualify = pd.get('qualify',''),
                               directed = False, force_qualify = True )
        else :
            temporaryds.parse( pdsf, fd = I.fd_biogrid, qualify = pd.get('qualify',''),
                               directed = False, force_qualify = True )

        print( 'saving public dataset' + pd['infilename'])
        sio = StringIO()
        temporaryds.save( sio, edges = { e for e in temporaryds.edges.values() 
                                         if e.weight >= pd.get('minweight',0) and e.totalscore >= pd.get('minscore',0) });
        sio.seek(0)
        print( 'imiport public data into network ' )
        nwdata.load_from( sio )
        sio.close()
        pdsf.close()

def networkwideRescue( edgeset, c ):
    # network-wide degree rescue
    # this block puts all edges with appropriate quals in the 
    # valid_qual_edges set, which will be used to calculate degrees
    # if NWD is true. Otherwise this set needs to be recreated
    # depending on which baits a protein binds

    if c['valid_quals'] == 'default' : 
        vqe = { e.key for e in edgeset if e.qual not in {'bg','bp'} }
    elif c['valid_quals'] == 'all' : 
        vqe = { e.key for e in edgeset }
    elif type(c['valid_quals']) is set : 
        vqe = { e.key for e in edgeset if e.qual in c['valid_quals'] }; 
    else : 
        vqe = set()
        
    return vqe


def rescueListNodes( c ):
    # read the rescue file
    if c[ 'rescue_f' ] :
        with open( c['rescue_f'] ) as resc:
            rescued = set( resc.read().splitlines()) 
    else :
        rescued = {}

    return rescued

def rescueEdgesByPublic( nwdata, c ):
    
    # public edges used to evaluate binding to baits

    # pass 1 edges include if ...:
    # this filtration will now be taken care of above, immediately after 
    #parsing
    
    reinforcing_edges = set()
    edges_pass1       = set()
    for e in list( nwdata.edges.values()) : 
        if e.key in c['zfhits_joint'] : 
            # passed filtering
            edges_pass1.add(e)

        elif {e.to.key,e.whence.key}.issubset( c['node_pass1_all'] ) and e.qual in\
             { pd.get('qualify','') for pd in c['public_dicts'] } and e.to.key != e.whence.key : 
            edges_pass1.add(e)
            reinforcing_edges.add(e)

    c['reinforcing_edges'] = reinforcing_edges
    c['edges_pass1']       = edges_pass1

    
def secondaryFiltration( nwdata, c ):

    rescueEdgesByPublic( nwdata, c )
    vqe            = networkwideRescue( nwdata.edges.values(), c ) if c['nwd'] else set( )
    rescued        = rescueListNodes( c )
    nnodes_rescued = 0
    nodes_pass2    = set()    

    # vet every node in the network
    for nk in c['node_pass1_all'] :
        if nk in c['node_pass1_strong'] :
            # keep node if it passed strong filter
            nodes_pass2.add(nk)
        elif b4us(nk) in rescued or afterus(nk) in rescued :
            # keep the node if it is on the rescue list:
            nodes_pass2.add(nk)
            nnodes_rescued  += 1
        elif any([ nwdata.nodes[nk].binds( bk, within_edge_set = c['reinforcing_edges']) and \
                   nwdata.nodes[nk].binds( bk, within_edge_set = c['zfhits_joint'] ) 
                   for bk in c['baitkeys'] ]) : 
            # if the node binds a bait (bk) with both a publicly recognized edge and
            # an observed weak edge, make sure that node stays
            nodes_pass2.add(nk)
            nnodes_rescued  += 1            
        elif not c['nwd'] :
            # rationale for this filter:
            edges_this_node = { e.key for es in nwdata.nodes[nk].edges.values() for e in es }\
                & ( c['edges_pass1'] | c['reinforcing_edges'] ); 
            # possibly ok edges
            vqe = networkwideRescue( edges_this_node, c )
            # should this be renamed NON-network wide rescue?

            for bk in c['baitkeys'] : 
                partners_this_node={ n for n in nwdata.nodes[nk].partners.values() if\
                                     nwdata.nodes[nk].binds( n, within_edge_set = vqe) and \
                                     nwdata.nodes[bk].binds( n, within_edge_set = c['edges_pass1'] ) }

                if len(partners_this_node) >= c['rescue_deg'] : 
                    nodes_pass2.add(nk)
                    nnodes_rescued  += 1                                
                    break

        elif c['nwd'] and nwdata.nodes[nk].nneighbors(within_edge_set = vqe) >= c['rescue_deg'] :
            # rationale for this filter:
            nodes_pass2.add(nk)
            nnodes_rescued  += 1            

    edges_pass2 = set()            
    for bk in c['baitkeys'] : 

        real_partners_this_node = { nk for nk in nodes_pass2 \
                                    if nwdata.nodes[bk].binds( nk, within_edge_set = c['edges_pass1'] | c['reinforcing_edges']) }

        real_partners_this_node.add(bk)

        edges_pass2    |= { e for e in c['edges_pass1'] | c['reinforcing_edges']
                            if {e.to.key,e.whence.key}.issubset(real_partners_this_node) }

    c['nodes_pass2']    = nodes_pass2
    c['edges_pass2']    = edges_pass2
    c['nnodes_rescued'] = nnodes_rescued
    

unitransform = lambda x : 7.0 if x==0.0 else -1 * np.log10(x)

def makeOutput( nwdata, c ):
        
    ie.print_springs( c['edges_pass2'], print_headers = True, print_weights = True, transform_scores = unitransform,
                      print_quals = True, fname = c['outfilename'] ,print_pps=True)

    sys.stdout.write("Nodes:\n  Pass 1: {}\n    Strong: {}\n  Pass 2: {}\n    Rescued : {}\n\n".
                     format( len( c['node_pass1_all'] ), len( c['node_pass1_strong'] ), len( c['nodes_pass2'] ), c['nnodes_rescued'] ))

    edgequals_1 = { e.qual for e in c['edges_pass1'] }
    edgequals_2 = { e.qual for e in c['edges_pass2'] }

    sys.stdout.write("Edges:\n  Pass 1: {}\n".format( len( c['edges_pass1'] )))
    for eq in edgequals_1 : 
        sys.stdout.write("{: <12}: {}\n".format( eq, len({ e for e in c['edges_pass1'] if e.qual == eq })))

    sys.stdout.write("\n  Pass 2: {}\n".format( len( c['edges_pass2'] )))
    for eq in edgequals_2 : 
        sys.stdout.write("{: <12}: {}\n".format( eq, len({ e for e in c['edges_pass2'] if e.qual == eq })))

    if 'idbfilename' in c :
        nwdata.save( c['idbfilename'], nodes = c['nodes_pass2'], edges = c['edges_pass2'] )


def createNetwork( yamlfile ) :

    readYAMLfile( yamlfile, config )
    loadObjects( config )
    
    theds = I.dataSet(n_filter = config['node_filter'], debug = DB)

    readInDatasets( theds, config )

    # filter experimental data by background dists 
    filterNodesByBackground( theds, config )

    readPublicDatasets( theds, config )

    secondaryFiltration( theds, config )

    makeOutput( theds, config )

    return theds,config
    
if __name__ == "__main__":

    theds,config=createNetwork( sys.argv[1] )

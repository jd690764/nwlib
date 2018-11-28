#!/usr/env/python -c
import argparse
import sys
import os
from os import listdir
from io import StringIO
import colorama
import numpy as np
import yaml
import re
from subprocess import call
import numpy
from lib import interactors as I
from lib import interactors_extras as ie
from lib import rbase
from lib.markutils import b4us,afterus
from lib import config as cf
from lib import utils
import getpass
from shutil import copyfile
import pprint
pp = pprint.PrettyPrinter(indent=4)

DB           = False
deprecation  = colorama.Fore.RED+" DEPRECATED 21 APRIL 2016 "+colorama.Fore.RESET
config       = { 'ifiles' : cf.ifilesPath,
                 'publicDatadir' : cf.publicDataPath,
                 'youtfilename' : '',
                 'outfilename'  : '',
                 'rescue_f'     : None }

scruser      = getpass.getuser()

def tokey_single(c, s) :
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

def tokey( c, s ):
    if type(s) is list:
        keys = list()
        for item in s:
            keys.append( tokey_single( c, item ))
        return keys
    else:
        return tokey_single( c, s )    

    
def loadObjects( c ):

    if c['organism'] == 'human' :
        rbase.load('hmg')
    elif c['organism'] == 'mouse' :
        rbase.load('mmg')

        
def readYAMLfile( yamlfile, c ) :

    with open( yamlfile ) as infile :
        yout          = yaml.load(infile.read())

    c['ds_dicts']     = yout['datasets']

    c['organism']      = yout['options'].get('organism','human')
    c['ALPHA_HI']      = yout['options'].get('alpha_hi',0.01)
    c['ALPHA_LO']      = yout['options'].get('alpha_lo',0.05)
    c['FLOOR_HI']      = yout['options'].get('floor_hi',4)
    c['FLOOR_LO']      = yout['options'].get('floor_lo',4)
    c['CORRL_HI']      = yout['options'].get('corrl_hi',0.70)
    c['CORRL_LO']      = yout['options'].get('corrl_lo',0.70)
    c['rescue_deg']    = yout['options'].get('degree_to_rescue',1)
    c['valid_quals']   = yout['options'].get('valid_degree_quals',{'wt',})
    c['mt_method']     = yout['options'].get('mt_method','fdr_bh')
    c['filter_by']     = yout['options'].get('filter_by','avgp')    
    c['node_filter']   = yout['options'].get('node_filter', None)
    c['iact_filter']   = yout['options'].get('iact_filter', None) 
    c['rescueNetWide'] = yout['options'].get('network-wide_degree',False)
    c['rescueLocal']   = yout['options'].get('rescue_local', False)
    c['rescueByComp']  = yout['options'].get('rescue_by_complex',False)
    c['rescueAll']     = yout['options'].get('rescue_all', False)
    c['public_dicts']  = yout['public']

    refFiles          = dict()
    for dsd in c['ds_dicts'] : 
        lost_files    = 0
        if not os.path.isfile('./'+dsd['infilename']) and \
           not os.path.isfile(c['ifiles']+dsd['infilename']) and \
           not os.path.isfile(dsd['infilename']): 
            print('File '+dsd['infilename']+' not found.')
            lost_files += 1
        if dsd['control'] not in refFiles:
            pFile     = ie.createReferenceFile( dsd['control'], overWrite = False)
            refFiles[dsd['control']] = pFile
            dsd['control_ori']       = dsd['control']            
            dsd['control']           = pFile
            
    if lost_files > 0 : 
        raise IOError

    # this should be a list of dicts
    # each dict (1/public file) could have the fields : infilename qualify convert misncore minweight bait [] radius 
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
        

def readInDatasets( nwdata, c, dict_to_use, keys ):

    # parse datasets
    baitkeys      = list()

    # iterating through the dataset items
    for dsd in c[ dict_to_use ] : 

        if dsd['infilename'] in os.listdir('.') : 
            dsf = open(dsd['infilename']) ; 
        elif os.path.isfile( c['ifiles'] + dsd['infilename'] ): 
            dsf = open(c['ifiles'] + dsd['infilename']) ; 
        elif os.path.isfile( dsd['infilename'] ) : 
            dsf = open(dsd['infilename']) ; 
        else :
            # problem
            continue

        convert  = dsd.get( 'convert', None )
        nwdata.parse( dsf, fd = I.fdms, convert = convert, qualify = dsd.get('qualify',''),
                      floor = c['FLOOR_HI'], directed = True, user = scruser )

        k = tokey(c, dsd['bait'])
        if type(k) is list: 
            baitkeys.extend( k )
        else:
            baitkeys.append( k )

        dsf.close()

    c[ keys ] = baitkeys
    

def filterNodesByBackground( nwdata, c ):

    zfhits_strong          = set( )
    zfhits_weak            = set( )
    #import pickle
    #sys.setrecursionlimit(1200)
    #pickle.dump(nwdata, open('/usr/local/share/py/djscripts/tmp/nwdata.pk', 'wb'))
    #pickle.dump(c, open('/usr/local/share/py/djscripts/tmp/nwdata_conf.pk', 'wb'))
    
    for dsd in c['ds_dicts'] :
        print(dsd['infilename'])
        zfhits_strong     |= ie.madfilter_corr( nwdata, dsd['control'], tokey(c, dsd['bait']),
                            convert = dsd.get('convert', None), 
                            qual = dsd.get('qualify'), directed = True, alpha = c['ALPHA_HI'],
                            maxcorr = c['CORRL_HI'], debug = True, method = c['mt_method'] )
        zfhits_weak       |= ie.madfilter_corr( nwdata, dsd['control'], tokey(c, dsd['bait']),
                            convert = dsd.get('convert', None),
                            qual = dsd.get('qualify'), directed = True, alpha = c['ALPHA_LO'],
                            maxcorr = c['CORRL_LO'], debug = DB, method = c['mt_method'] )

    store_first_pass_data( nwdata, c, zfhits_strong, zfhits_weak )
    
        
def store_first_pass_data( nwdata, c, strong_hits, weak_hits ):
    # expt bait to node edges, that are deemed significant
    joint_hits             = strong_hits | weak_hits

    # pass 1 nodes : every node at the end of one of the validated edges
    # remove nodes that don't have experimental edges pointing to them
    node_pass1_all         = { nk for ek in joint_hits  for nk in { nwdata.edges[ek].to.key, nwdata.edges[ek].whence.key}}
    node_pass1_strong      = { nk for ek in strong_hits for nk in { nwdata.edges[ek].to.key, nwdata.edges[ek].whence.key}}

    c['joint_hits']        = joint_hits
    c['node_pass1_all']    = node_pass1_all
    c['node_pass1_strong'] = node_pass1_strong

def trim( nw, nodes = [], rad = 1 ) :

    final_edge_set          = set()

    if len(nodes) > 0 and rad >= 1:
        base_node_set       = set( nodes )
        final_node_set      = base_node_set.copy()
        while rad > 0:
            final_node_set |= {e.whence.official for e in nw.edges.values() if e.to.official in base_node_set }
            final_node_set |= {e.to.official for e in nw.edges.values() if e.whence.official in base_node_set }
            rad             = rad -1
            base_node_set   = final_node_set.copy()

        #print( 'nodes: ' + '+'.join(nodes) + '\n')
        #print( 'final_node_set: ' +  '+'.join(map(str, list(final_node_set)[0:10])))
        #print( 'edge_values: ' + '+'.join(map(str, list({e.to.official for e in nw.edges.values()})[0:10])) + '\n')
        final_edge_set     = {e for e in nw.edges.values() if {e.to.official,e.whence.official}.issubset(final_node_set) }
        #print('edges in set: ' + str(len(final_edge_set)))
    return final_node_set
    
def readPublicDatasets( nwdata, c ):

    if type(c['public_dicts']) is not list or len(c['public_dicts']) == 0:
        return True

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
        convert     = pd.get('convert', None)
        
        temporaryds.parse( pdsf, fd = I.fd_biogrid, convert = convert, qualify = pd.get('qualify',''),
                               directed = False, force_qualify = True, user = scruser )

        sio = StringIO()

        print( 'saving public dataset' + pd['infilename'])
        # filter nodes of the dataset if bait is defined
        if 'bait' in pd:
            node_set = trim( temporaryds, pd.get('bait'), pd.get( 'radius', 1 ))
            temporaryds.save( sio, nodes = node_set );
        else :
            edge_set =  { e for e in temporaryds.edges.values() if e.weight >= pd.get('minweight',0) and e.totalscore >= pd.get('minscore',0) }
            temporaryds.save( sio, edges = edge_set );

        sio.seek(0)
        print( 'import public data into network ' )
        nwdata.load_from( sio, scruser )
        sio.close()
        pdsf.close()
        #temporaryds.save( cf.djPath + 'tmp/adcy3_tmp.txt', edges = edge_set );
        #print('number of nodes: ' + str(len([n for n in nwdata.nodes])))
        
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
        if e.key in c['joint_hits'] : 
            # passed filtering
            edges_pass1.add(e)

        elif {e.to.key,e.whence.key}.issubset( c['node_pass1_all'] )\
             and e.qual in { pd.get('qualify','') for pd in c['public_dicts'] }\
             and e.to.key != e.whence.key : 
            edges_pass1.add(e)
            reinforcing_edges.add(e)

    c['reinforcing_edges'] = reinforcing_edges
    c['edges_pass1']       = edges_pass1

def getComplexes( c ):
    # read in complexes
    complexes = dict()
    with open( cf.complexByGeneFile, 'rt' ) as cbg:
        for line in cbg:
            v = line.rstrip().split('\t')
            complexes[ v[0] ] = set( v[1].split('|') )
    return complexes
    
def secondaryFiltration( nwdata, c ):

    rescueEdgesByPublic( nwdata, c )
    vqe            = networkwideRescue( nwdata.edges.values(), c ) if c['rescueNetWide'] else set( )
    complexes      = getComplexes( c ) if c['rescueByComp'] else dict()
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
                   nwdata.nodes[nk].binds( bk, within_edge_set = c['joint_hits'] ) 
                   for bk in c['baitkeys'] ]) : 
            # keep the node if it binds a bait (bk) with both a public edge and
            # an experimental weak edge
            nodes_pass2.add(nk)
            nnodes_rescued  += 1            

        elif c['rescueByComp'] and len([ on for on in c['node_pass1_all']
                                         if nk != on
                                         and b4us(nk) in complexes
                                         and b4us(on) in complexes
                                         and bool( complexes[b4us(nk)] & complexes[b4us(on)])]) >= c['rescue_deg']:
            # keep node if there is another member of a complex among 'node_pass1_all'
            nodes_pass2.add( nk )
            nnodes_rescued += 1

        elif c['rescueLocal'] :
            # keep the node if it has at least 'rescue_deg' neighbors that themselves
            # bind a bait with edges within 'edges_pass1'
            # call this 'local rescue'
            edges_this_node = { e.key for es in nwdata.nodes[nk].edges.values() for e in es } & c['edges_pass1']
            # possibly ok edges
            vqe = networkwideRescue( edges_this_node, c )
            for bk in c['baitkeys'] : 
                partners_this_node={ n for n in nwdata.nodes[nk].partners.values() if\
                                     nwdata.nodes[nk].binds( n, within_edge_set = vqe) and \
                                     nwdata.nodes[bk].binds( n, within_edge_set = c['edges_pass1'] ) }

                if len(partners_this_node) >= c['rescue_deg'] : 
                    nodes_pass2.add(nk)
                    nnodes_rescued  += 1                                
                    break

        elif c['rescueNetWide'] and nwdata.nodes[nk].nneighbors(within_edge_set = vqe) >= c['rescue_deg'] :
            # keep node if it has at least 'rescue_deg' neighbors within all qual edges
            nodes_pass2.add(nk)
            nnodes_rescued += 1

    edges_pass2 = set()            
    for bk in c['baitkeys'] : 
        if bk in nwdata.nodes:
            real_partners_this_node = { nk for nk in nodes_pass2 \
                                        if nwdata.nodes[bk].binds( nk, within_edge_set = c['edges_pass1'] ) }
            real_partners_this_node.add(bk)

            edges_pass2    |= { e for e in c['edges_pass1'] if {e.to.key,e.whence.key}.issubset(real_partners_this_node) }
        else:
            print( bk + ' is not among the network nodes!' )
            
    c['nodes_pass2']    = nodes_pass2
    c['edges_pass2']    = edges_pass2
    c['nnodes_rescued'] = nnodes_rescued
    

unitransform = lambda x : 7.0 if x==0.0 else -1 * np.log10(x)

def makeOutput( nwdata, c ):

    ie.print_springs( c['edges_pass2'], print_headers = True, print_weights = True, transform_scores = unitransform,
                      print_quals = True, fname = c['outfilename'] ,print_pps=True, print_bkg = True )

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

def makeSaintInputFiles( nwdata, cdata, c, add_replicate = False ):

    c['interfile'] = re.sub(r'.cyt', r'.inter.txt', c['outfilename'])
    c['preysfile'] = re.sub(r'.cyt', r'.preys.txt', c['outfilename'])
    c['baitsfile'] = re.sub(r'.cyt', r'.baits.txt', c['outfilename'])
    c['listfile']  = re.sub(r'.cyt', r'.list.txt', c['outfilename'])
    
    inter = dict() # all the interactions to be put in output file
    preys = dict() # all preys in non-control datasets
    baits = dict() # all baits in samples/controls
    numpy.random.seed( 123 )
        
    for dsd in c['ds_dicts'] :

        sample = re.sub(r'.i', r'', dsd['infilename'])
        # the bait_qualify value is going to be the 'bait' 
        baits[ dsd['infilename']] = dsd.get('bait') + '_' + dsd.get('qualify') + "\t" + 'T'
        if add_replicate:
            baits[ dsd['infilename'] + '.r' ] = dsd.get('bait') + '_' + dsd.get('qualify') + "\t" + 'T'
        bkey = tokey(c, dsd['bait'])
        for e in { e for es in nwdata.nodes[ bkey ].edges.values() for e in es  } :
            if dsd.get('qualify') and dsd.get('qualify') != e.qual:
                continue
            elif e.to.key == bkey or e.to.key == 'PSEUDO_00':
                continue

            prey = re.sub(r'^(.+)_\d+$', '\\1', e.to.key)
            #sc   = int(sum([ int(re.sub(r'^.+_(\d+)$', '\\1', i.tags)) for i in e.interactions ]) / len(e.interactions))
            sc   = [ int(re.sub(r'^.+_(\d+)$', '\\1', i.tags)) for i in e.interactions if sample in i.interID ]
            aal  = sum([ float(re.sub( r'^.+len_(.+)_raw.*$', '\\1', i.tags )) for i in e.interactions ]) / len(e.interactions)
            i    = dsd['infilename'] + "\t" + dsd.get('bait') + '_' + e.qual + "\t" + prey 
            ir   = dsd['infilename'] + '.r'  + "\t" + dsd.get('bait') + '_' + e.qual + "\t" + prey

            #print('+'.join(map(str,sc)) + ' ' + str(type(sc)), sample, str(len(sc)))
            if len(sc) > 0 and sc[0] > 1: # skip preys with spectral count = 1
                preys[ prey ]   = aal
                if i not in inter:
                    inter[ i ]  = list()
                if add_replicate and ir not in inter:
                    inter[ ir ] = list()
                inter[ i ].append( sc[0] )
                if add_replicate:
                    scr = int(numpy.random.normal( sc , 10, size = None))
                    scr = scr if scr > 0 else 0
                    inter[ ir ].append( scr )

    counter   = 1
    for dsd in c['c_dicts'] :
        ckey          = 'C_' + str( counter )
        counter       = counter + 1
        baits[ ckey ] = ckey + "\t" + 'C'
        bkey          = tokey(c, dsd['bait'])
        for e in { e for es in cdata.nodes[ bkey ].edges.values() for e in es  } :
            if dsd.get('qualify') and dsd.get('qualify') != e.qual:
                continue
            elif e.to.key == bkey or e.to.key == 'PSEUDO_00':
                continue

            prey = re.sub(r'^(.+)_\d+$', '\\1', e.to.key)
            sc   = int(sum([ int(re.sub(r'^.+_(\d+)$', '\\1', i.tags)) for i in e.interactions ]) / len(e.interactions))
            i    = ckey + "\t" + ckey + "\t" + prey 

            if prey in preys:
                if i not in inter:
                    inter[ i ] = list()
                inter[ i].append( sc )

        for p in preys:
            ip   = ckey + "\t" + ckey + "\t" + p 
            if ip not in inter:
                inter[ ip ] = [ 0 ]
    
    with open( c['interfile'], 'wt') as oh:
        for k, v in sorted(inter.items()):
            oh.write( k + '\t' + ','.join(map(str, v)) + "\n")
    with open( c['preysfile'], 'wt') as oh:
        for k, v in sorted(preys.items()):
            oh.write( k + '\t' + str( v ) + "\n")
    with open( c['baitsfile'], 'wt') as oh:
        for k, v in sorted(baits.items()):
            oh.write( k + '\t' + v + "\n")
            
def readControls( c ):

    cntrl = I.dataSet(n_filter = config['node_filter'], debug = DB)
    for dsd in c['ds_dicts'] :
        cfile = dsd['control_ori']
        # assumed that this is the same in all datasets and
        # it is a file of a list of ifilenames
        break

    c_dicts = list()
    with open( cfile, 'rt' ) as fh:
        for line in fh:
            ifname = line.rstrip()
            bait  = line.split('_')[0]
            c_dicts.append( { 'infilename': ifname, 'bait': bait, 'qualify': ifname })

    c['c_dicts'] = c_dicts
    readInDatasets( cntrl, c, 'c_dicts', 'cntrlkeys' )
    return cntrl

def filterSaintData( ds, c, saint, cutoff, cutoff_value, update = False ):

    selected_edges = set( )
    
    # assign saint scores to appropriate edges   
    for dsd in c['ds_dicts']:
        bait = dsd['bait']
        ek_prey, x, y = ie.dataset_edges_for_bait( ds, baitkey = tokey(c, bait), qual = dsd.get('qualify', ''), directed = True, debug = DB )
        for ek in ek_prey.keys() :
            values = [ 0, 1 ]
            e      = ds.edges[ek]
            qual   = e.qual
            prey   = e.to.official
            # assign FDR to p and avgP to a new attr, avgP
            if bait in saint and prey in saint[bait] and qual in saint[bait][prey] and isinstance( saint[bait][prey][qual], list ):
                values = saint[bait][prey][qual]
                if update:
                    ds.edges[ek].p    = saint[bait][prey][qual][ 1 ]
                    ds.edges[ek].avgP = saint[bait][prey][qual][ 0 ]
                    ds.edges[ek].bkg  = saint[bait][prey][qual][ 2 ]
                    ds.edges[ek].spc  = saint[bait][prey][qual][ 3 ]
                    
            if cutoff == 'fdr' and values[1] < cutoff_value:
                selected_edges.update( [ek] )
            elif cutoff == 'avgp' and values[0] > cutoff_value:
                selected_edges.update( [ek] )

    return selected_edges
                
def scoreBySaintx( ds, c ):
    os.chdir('/usr/local/share/py/djscripts/tmp/')

    # network for control experiments
    cntrl = readControls( c )
    # create input files for saintExpress
    makeSaintInputFiles( ds, cntrl, c, add_replicate = False )
    # run saintExpress
    from subprocess import Popen, PIPE

    p = Popen( [ './../bin/SAINTexpress-spc ' + c['interfile'] + ' ' +  c['preysfile'] + ' ' +c['baitsfile'] ],
               cwd = '/usr/local/share/py/djscripts/tmp/',
               stdin=PIPE, stdout=PIPE, stderr=PIPE, shell = True )
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    rc = p.returncode
    if( err ):
        print( 'rc=', str(rc), "\n", output, "\n", err )    

    # the output is in the root directory
    copyfile( 'list.txt', c['listfile'] )
    # read in saint scores
    saint_data = dict()
    with open( c['listfile'] ) as lf:
        for line in lf:
            if re.search( r'^Bait', line ):
                continue
            fields = line.rstrip().split( '\t' )
            y      = fields[ 0 ].split( '_' )
            # save bait, prey, qual, avgP, FDR, ctrlCounts
            if y[0] not in saint_data:
                saint_data[ y[0] ] = dict()
            if fields[1] not in saint_data[ y[0] ]:
                saint_data[ y[0] ][ fields[1] ] = dict()
            if y[1] not in saint_data[ y[0] ][ fields[1] ]:
                saint_data[ y[0] ][ fields[1] ][ y[1] ] = list()
            saint_data[ y[0] ][ fields[1] ][ y[1] ].extend( [ float(fields[8]), float(fields[15]), fields[7], fields[3] ] )

    hits_strong  = filterSaintData( ds, c, saint_data, c['filter_by'], c['ALPHA_HI'], True )
    hits_weak    = filterSaintData( ds, c, saint_data, c['filter_by'], c['ALPHA_LO'], False )

    store_first_pass_data( ds, c, hits_strong, hits_weak )    

def makeGctFiles( cf ):

    ifiles = list()
    for ds in cf['ds_dicts']:
        ifiles.append( ds['infilename'] )

    utils.makeGct( ifiles, cf['outfilename'] + '.gct', scoreCol = 9 )

    
def createNetwork( yamlfile ) :

    readYAMLfile( yamlfile, config )
    loadObjects( config )

    theds = I.dataSet(n_filter = config['node_filter'], debug = DB)

    readInDatasets( theds, config, 'ds_dicts', 'baitkeys' )
    
    if config['mt_method'] == 'fdr_bh':
        # filter experimental data by background dists 
        filterNodesByBackground( theds, config )
        
    elif config['mt_method'] == 'saintx':
        scoreBySaintx( theds, config )
    else:
        config['joint_hits']        = set(theds.edges.keys())
        config['node_pass1_all']    = set(theds.nodes.keys())
        config['node_pass1_strong'] = set(theds.nodes.keys())

    readPublicDatasets( theds, config )

    if not config['rescueAll']:
        secondaryFiltration( theds, config )

    else :
        
        config['nodes_pass2']    = config['node_pass1_all']
        config['edges_pass1']    = set()

        for e in list( theds.edges.values()) : 
            if {e.to.key,e.whence.key}.issubset( config['node_pass1_all'] ) and e.to.key != e.whence.key : 
                config['edges_pass1'].add(e)
        
        config['edges_pass2']    = config['edges_pass1']
        config['nnodes_rescued'] = 0

    
    makeOutput( theds, config )
    
    return theds,config
    
if __name__ == "__main__":

    theds,config=createNetwork( sys.argv[1] )

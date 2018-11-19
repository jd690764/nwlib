#foo
import sys
import re
from urllib.error import HTTPError
import numpy as np
from os.path import isfile
from warnings import warn
import pickle
import os
from scipy.stats.distributions import norm, truncnorm
from statsmodels.sandbox.stats.multicomp import multipletests
from mpmath import mp as mpmns

from lib import interactors as I 
from lib import markutils as mu
from lib import rbase
from lib import config as cf

sow=sys.stdout.write
sew=sys.stderr.write

CONTROL_FILES  = cf.controlFiles
PSEUDO_DEFAULT = 1e-5

def topgroup(keykeys,wns,wes=None,remove_empty_groups=True,print_summary=True) : 
    """
    interactors_extras.topgroup(keykeys,wns) 

    searches the dataset for the nodes whose keys are listed in keykeys, then groups
    all remaining nodes in the dataset by which nodes they have edges to, either using
    all edges related to those nodes or strictly those in wes.

    beginning 23 april 2014, items in keykeys can be sets of genes 

    While not strictly necessary, it is STRONGLY recommended that wns be an ordered iterable such
    that the tuples in dict 1 (see below) make sense to you.

    Returns two dicts.
    Dict 1 : keyed by  a tuple of bools indicating whether the nodes
        in that topology group bind to the node represented by that place in the tuple or no.

        I.e. if you submitted wns=  { node_a,node_b,node_c }, then topology group 0, where the nodes bind
            none of the query nodes, should be { (False,False,False) : 0 }

        Values are the number of the topology group, and the key to dict 2

    Dict 2 : keyed by node key,values are also topology group.


    Be careful with keykeys-- the number of possible topology groups is 2**(len(keykeys))
    """

    ####### create topology groups -- list of lists
    # recursively creates all possible lists of bools with a given length
    warn(DeprecationWarning('I (Mark) can\'t vouch for topgroup right now 3/1/2016') )
    def rboollist(masterlist,inlist,length) : 

        wt=inlist + [True ];
        wf=inlist + [False] ; 

        if len(wt) == length  : 
            masterlist.extend([wt,wf]) ; 

        else : 
            rboollist(masterlist,wt,length) ; 
            rboollist(masterlist,wf,length) ; 

    masterlist=list() ; 
    inlist=list() ; 

    rboollist(masterlist,inlist,len(keykeys)) ; 

    uberlist=list() ; 

    d1=dict() ;
    d2=dict() ;

    for l in masterlist :

        uberlist.append(tuple(l)) ; 

    del masterlist ; 

    # d1 : key is tuple of bools, value is group id # 
    # rd1 : key is group id #, value is tuple of bools 
    d1=dict(list(zip(uberlist,list(range(0,len(uberlist)))))) ; 
    rd1=dict(list(zip(list(range(0,len(uberlist))),uberlist))) ; 

    # d2 : initialized above : key is node key (from wns), value is topology tuple

    for n in wns : 

        boollist=list() ; 

        for k in keykeys : 

            if type(k) is str : 
                if wes :
                    if n.edges.get(k) in wes : 
                        boollist.append(True) ; 
                    else:
                        boollist.append(False) ; 
                elif k in n.edges : 
                    boollist.append(True) ; 
                else : 
                    boollist.append(False) ; 
            elif type(k) is set : 
                if wes : 
                    for subk in k : 
                        if n.edges.get(subk) in wes : 
                            boollist.append(True) ; 
                            break ; 
                    else : 
                        boollist.append(False) ; 
                else : 
                    for subk in k : 
                        if n.edges.get(subk) : 
                            boollist.append(True) ; 
                            break ; 
                    else : 
                        boollist.append(False) ; 

        booltup=tuple(boollist) ; 

        d2.update({ n.key : d1[booltup] }) ; 

    popgroups=set() ; 
    sd2v=set(d2.values())

    for k,v in list(d1.items()) : 
        # search through d1 ( topology tuple : # ) 
        if v in sd2v : 
            # if the group # is represented in d2, add it to 'populated groups'
            popgroups.add(v) ;

    if remove_empty_groups : 

        for k,v in list(d1.items()) : 
            if v not in popgroups : 
                d1.pop(k) ; 

    if print_summary : 

        sow('SUMMARY : \n') ; 
        sow('    Possible groups : {} \n'.format(2**len(keykeys))) ; 
        sow('    Populated groups : {} \n'.format(len(popgroups))) ; 

        sow('    GROUP:\tMEMBERS:\tTOPOLOGY:\n')

        popdict=dict() ;

        for v in popgroups :
            popdict.update({ v : len([x for x in list(d2.keys()) if d2.get(x) == v])}) ; 

        for p in popdict : 
            sow('    {}\t{}\t{}\n'.format(p,popdict[p],rd1[p])) ; 

    return (d1,d2) ; 


def prep_preppi(preppi_fobj,outfile,h2m=False,completionist=False) : 
    warn(DeprecationWarning('I (Mark) can\'t vouch for prep_preppi right now 3/1/2016') )

    #import rbase
    rbase.load('dup') ;
    rbase.load('hmg') ; 
    from lib.markutils import date6,getLineCount,waitbar
    from peputils import bestPep as bp
    from engarde import fetchSw as sw

    # creates an interactors-style object from downloaded PrePPI table

    idprefix='PrePPI_'+date6()  ;

    if h2m : 
        rbase.load('h2m') ; 
        theorg='10090'
    else : 
        theorg='9606'

    dup=rbase.dup ;
    h2m=rbase.h2m ; 
    hmg=rbase.hmg ; 

    swaccmap=dict() ; 

    nlines=getLineCount(preppi_fobj) ; 

    i=0 ; 
    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(\
     'ID','LEFT','RIGHT','LR','ORG','ORG','EID_LEFT','EID_RIGHT','ID_PREFIX','')) ; 

    for line in preppi_fobj : 

        line_failed=False
        waitbar(80 * i / nlines,80,showPct=True) ; 

        linel=line.strip().split('\t') ; 

        if linel[10] != 'NA' : 
            i += 1; 
            continue ; 

        eids=[ '','' ];

        for n in [0,1] :
            eids[n]= swaccmap[linel[n]] if linel[n] in swaccmap else \
                  dup[linel[n]] if linel[n] in dup else \
                  hmg['SwissProt'][linel[n]]['EID'] if type(hmg['SwissProt'].get(linel[n])) is dict else \
                  hmg['SwissProt'][linel[n]][0]['EID'] if type(hmg['SwissProt'].get(linel[n])) is list else \
                  None ; 

            if not eids[n] and completionist : 
                try : 
                    swacc,seq=sw(linel[n],asTuple=True) ; 
                    pepacc=bp(seq) ; 
                    if pepacc not in hmg['Peptide'] : 
                        #print 'Failed swissprot-to-entrez conversion.'
                        line_failed=True ; 
                        break ; 
                except (HTTPError,IOError) :
                    line_failed=True ; 
                    break ; 
                eids[n]=hmg['Peptide'][pepacc]['EID'] ; 
            elif not eids[n] : 
                line_failed=True ; 
                break ; 

            if h2m : 
                eids[n] =   h2m[eids[n]] if type(h2m[eids[n]]) is str else \
                            h2m[eids[n]][0] if type(h2m[eids[n]]) is list else \
                            None ; 
                if not eids[n] : 
                    #print 'Failed human-to-mouse conversion.'
                    line_failed=True ; 
                    break ; 

            if eids[n] and linel[n] not in swaccmap : 
                swaccmap.update({ linel[n] : eids[n] }) ; 

        if line_failed : 
            i+= 1 ; 
            continue ; 

        syms=[ hmg['EID'][x]['Symbol'] for x in eids ];

        outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(\
         idprefix+'_'+repr(i),syms[0],syms[1],linel[9],theorg,theorg,eids[0],eids[1],idprefix,'')) ; 
        # added blank string @ end for the tags field

        i += 1;

def bgsubtract_outkeys(query_dses,bg_dses,query_baits,bg_baits,cutoff=0.0,as_dict=False,debug=False) : 
    """
        bgsubtract(query_dses,bg_dsses,cutoff,metric='rel_dif')  ;

        returns tuple of lists (or dicts) of nodes in query_dses for which the
        metric for those nodes >= cutoff. 

        If the argument as_dict is provided, cutoff is ignored and all keys
        are returned as keys to a dict where the values are the relative
        difference of the query signal and the bait signal (over the bait signal) 

        The metric is presumed to incorporate the signal of the corresponding nodes
        in the bg_dsses.

        Cutoff can be an integer, float or ordered container. If a tuple, it must have the same
        length as bg_dses and correspond to different cutoffs for different datasets.
    """
    warn(DeprecationWarning('I (Mark) can\'t vouch for bgsubtract_outkeys right now 3/1/2016') )

    try :
        next(iter(cutoff)) ;
    except TypeError : 
        cutoff=[cutoff] ; 

    try :
        next(iter(query_dses)) ;
    except TypeError : 
        query_dses=[query_dses] ; 
    
    try :
        next(iter(bg_dses)) ;
    except TypeError : 
        bg_dses=[bg_dses] ; 

    if type(query_baits) is str : 
        query_baits=[query_baits] ; 

    if type(bg_baits) is str : 
        bg_baits=[bg_baits] ; 

    #if len(cutoff) > 1 and len(cutoff) != len(bg_dses) : 
    #    sys.stderr.write('ERROR : If length of cutoff >=1, its length must match # of background dses.\n') ; 
    #    raise IndexError ; 

    #if len(query_dses) != len(query_baits) and len(query_dses) > 1 :
    #    sys.stderr.write('ERROR : Length of query_dses > 1 and query_baits must match.\n') ; 
    #    raise IndexError ; 
    #elif len(query_dses) != len(query_baits) : 
    #    # this is hackish in the extreme--
    #    # if you have ONE dataset (like a merged MS set) with which you want to interrogate multiple baits,
    #    # we make a list whose length matches that of the baits where each element points to the same data set
    #    tmp=query_dses[0] ; 
    #    query_dses=list() ; 
    #    while len(query_dses ) < len(query_baits) : 
    #        query_dses.append(tmp) ; 
    #if len(query_dses) != len(query_baits) :
    #    sys.stderr.write('ERROR : Length of bg_dses and bg_baits must match.\n') ; 
    #    raise IndexError ; 

    #sync lengths of query, query_baits
    if len(query_dses) != len(query_baits) : 
        if debug: 
            sys.stderr.write('Syncing query_dses (length {}) and query_baits (length {}).\n'\
             .format(len(query_dses),len(query_baits)))
        if len(query_dses) == 1 : 
            query_dses=[ query_dses[0] for i in mu.indices(query_baits) ] ; 
        elif len(query_baits) ==  1 : 
            query_baits=[ query_baits[0] for i in mu.indices(query_dses) ] ; 
        else : 
            sys.stderr.write('Can\'t sync lengths of query lists where both have len > 1.\n') ; 
            raise IndexError ; 

    #sync lengths of bg_dses, query_baits
    if len(bg_dses) != len(bg_baits) : 
        if debug: 
            sys.stderr.write('Syncing bg_dses (length {}) and bg_baits (length {}).\n'\
             .format(len(bg_dses),len(bg_baits)))
        if len(bg_dses) == 1 : 
            bg_dses=[ bg_dses[0] for i in mu.indices(bg_baits) ] ; 
        elif len(bg_baits) ==  1 : 
            bg_baits=[ bg_baits[0] for i in mu.indices(bg_dses) ] ; 
        else : 
            sys.stderr.write('Can\'t sync lengths of background lists where both have len > 1.\n') ; 
            raise IndexError ; 

    if len(query_dses) != len(cutoff) : 
        if debug: 
            sys.stderr.write('Syncing cutoffs (length {}) and controls (length {}).\n'\
             .format(len(bg_dses),len(bg_baits)))
        if len(cutoff ) == 1 : 
            cutoff=[ cutoff[0] for i in mu.indices(bg_baits) ] ; 
        else : 
            sys.stderr.write('Can\'t sync length of cutoffs to bg_dses where both have len > 1.\n') ; 
            raise IndexError ; 


    pseudos=[ np.mean([e.meanscore for e in bgd.nodes[bgb].edges['PSEUDO_00'] ]) for bgd,bgb in zip(bg_dses,bg_baits) ]  ;

    print([ query_dses[i].infilenames for i in range(0,len(query_dses)) ]) ; 

    if debug : 
        sew('PSEUDOCOUNTS :') ;  
        for ps in pseudos : 
            sew('\t{:6.2}'.format(ps)) ; 
        sew('\n') ; 
        

    out_pretup=list() ; 

    for qd,qb in zip(query_dses,query_baits) : 

        if as_dict : 
            outdict=dict() ; 
        else : 
            outset=set() ; 

        for k in list(qd.nodes.keys()) : 

            if debug : 
                sew('NODE : {} in {} against bait {}\n'.format(k,query_dses.index(qd),qb)) ; 

            #!!
           #bg_sigs=[ bgd.nodes[bgb].edges[k].meanscore if bgd.nodes[bgb].edges.has_key(k) else bps \
           #           for bgd,bgb,bps in zip(bg_dses,bg_baits,pseudos) ]
            bg_sigs=list() ; 
            for bgd,bgb,bps in zip(bg_dses,bg_baits,pseudos) : 
                if k in bgd.nodes[bgb].edges : 
                    bg_sigs.append(max([ e.meanscore for e in bgd.nodes[bgb].edges[k] ])) ; 
                else:
                    bg_sigs.append(bps) ; 

            #!!
            #query_sig=qd.nodes[qb].edges[k].meanscore if qd.nodes[qb].edges.has_key(k) else 0.0 ; 
            query_sig=0.0 if k not in qd.nodes[qb].edges \
                      else max([ e.meanscore for e in qd.nodes[qb].edges[k]]) ;
                            

            if debug : 
                sew('QUERY SIG : {:6.2} \n'.format(query_sig)) ; 
                sew('BG SIGS : ') ; 
                for s in bg_sigs : 
                    sew('\t{:6.2}'.format(s)) ; 
                sew('\n') ; 

            rd=[ ( query_sig - bgs ) / bgs for bgs in bg_sigs ]  ;
            if debug : 
                sew('REL DIFFS: ') ; 
                for s in rd : 
                    sew('\t{:6.2}'.format(s)) ; 
                sew('\n') ; 
            

            if as_dict : 
                outdict.update({ k : rd }) ; 

            elif len(cutoff) == 1 :  
                # same cutoff for all bg sets
                for rdi in rd : 
                    if rdi < cutoff[0] : 
                        if debug : 
                            sew('    RD {:6.2} fails against cutoff {:6.2}\n'.format(rdi,cutoff[0])) ; 
                        break ; 
                else : 
                    if debug : 
                        sew('    Adding {}\n'.format(k)) ; 
                    outset.add(k) ;
            else : 
                # different cutoff for each bg set
                for rdi,cut in zip(rd,cutoff) : 
                    if rdi < cut : 
                        if debug : 
                            sew('    RD {:6.2} fails against cutoff {:6.2}\n'.format(rdi,cut)) ; 
                        break ; 
                else : 
                    if debug : 
                        sew('    Adding {}\n'.format(k)) ; 
                    outset.add(k) ; 

        if as_dict : 
            out_pretup.append( outdict ) ;
        else : 
            out_pretup.append( outset ) ;

    return tuple(out_pretup) ; 


def domenrich(query_ds,bg_dses,query_baits,bg_baits,taxon='9606',directed=None,qual=None,debug=False) : 
    warn(DeprecationWarning('I (Mark) can\'t vouch for domenrich right now 3/1/2016') )
    from numpy import log10,mean,std
    from lib.markutils import b4us,afterus

    if taxon == '9606' : 
        rbase.load('hsg') ; 
        rbase.load('hsp') ; 
        pept=rbase.hsp
        gene=rbase.hsg
    elif taxon == '10090' : 
        rbase.load('mmg') ; 
        rbase.load('mmp') ; 
        pept=rbase.mmp
        gene=rbase.mmg
    else : 
        sys.stderr.write('No other organisms ready yet.') ; 

    rbase.load('cdd') ; 

    domsigs =   dict() ; 
    dombg_means  =   dict() ;
    dombg_vars   =   dict() ; 

    counted_edges=set() ; 

    # get signals
    for bk in query_baits : 

        for k in list(query_ds.nodes[bk].partners.keys()) : 

            if afterus(k) not in gene['EID'] or gene['EID'][afterus(k)]['CDD'] is None : 
                continue ; 

            cdds={ rbase.cdd['ByPssm'][c]['Root']['Pssm']  \
             for c in gene['EID'][afterus(k)]['CDD'] if c in rbase.cdd['ByPssm']  } ; 

            for c in cdds : 

                for e in query_ds.nodes[bk].edges[k] : 
                    if (qual is None or e.qual in qual ) and \
                       (directed is None or (directed == e.directed and e.whence.key == bk)) and\
                        e.key not in counted_edges :

                        if c not in domsigs : 
                            domsigs.update({ c : 0.0 }) ; 

                        domsigs[c] += e.meanscore
                        counted_edges.add(e.key) ; 


    dombgs=dict()  
    # get backgrounds
    for bgbk,bgds in zip(bg_baits,bg_dses) : 
        for k in list(bgds.nodes[bgbk].partners.keys()) : 

            if afterus(k) not in gene['EID'] or gene['EID'][afterus(k)]['CDD'] is None : 
                continue ; 

            cdds={ rbase.cdd['ByPssm'][c]['Root']['Pssm']  \
             for c in gene['EID'][afterus(k)]['CDD'] if c in rbase.cdd['ByPssm']  } ; 

            for c in cdds : 
                if c not in dombgs and c in domsigs : 
                    dombgs.update({ c : [ sum([ e.meanscore for e in bgds.nodes[bgbk].edges[k] ]) ] }) ; 
                elif c in domsigs : 
                    dombgs[c].append(sum([ e.meanscore for e in bgds.nodes[bgbk].edges[k] ])) ; 
    
    # sum of scores for all pseudocounts in null set
    pseudos = [ sum([ e.meanscore for e in bgds.nodes[bgbk].edges['PSEUDO_00'] ]) 
                for bgbk,bgds in zip(bg_baits,bg_dses) ] ; 

    domenrich=dict() ; 
    for c in domsigs : 

        vals=dombgs.get(c,pseudos) ; 
        for i in range(0,len(bg_dses)-len(vals)) : 
            vals.append(0.0) ; 

        cmean=mean( vals ) ; 
        cstdv=std( vals ) ; 

        domenrich.update({ c : ( len(bg_baits) *  domsigs[c]  - len(query_baits) * cmean ) / cstdv }) ; 

    return domenrich


def bgsubtract(query_dses,bg_dses,query_baits,bg_baits,cutoff=0.0,as_dict=False,debug=False,selfs=False) : 
    """
        bgsubtract(query_dses,bg_dsses,cutoff,metric='rel_dif')  ;

        returns tuple of lists (or dicts) of EDGES in query_dses for which the
        metric for the non-bait nodes to the designated bait >= cutoff. 

        If the argument as_dict is provided, cutoff is ignored and all keys
        are returned as keys to a dict where the values are the relative
        difference of the query signal and the bait signal (over the bait signal) 

        The metric is presumed to incorporate the signal of the corresponding nodes
        in the bg_dsses.

        Cutoff can be an integer, float or ordered container. If a tuple, it must have the same
        length as bg_dses and correspond to different cutoffs for different datasets.
    """
    warn(DeprecationWarning('I (Mark) can\'t vouch for bgsubtract right now 3/1/2016') )

    try :
        next(iter(cutoff)) ;
    except TypeError : 
        cutoff=[cutoff] ; 

    try :
        next(iter(query_dses)) ;
    except TypeError : 
        query_dses=[query_dses] ; 
    
    try :
        next(iter(bg_dses)) ;
    except TypeError : 
        bg_dses=[bg_dses] ; 

    if type(query_baits) is str : 
        query_baits=[query_baits] ; 

    if type(bg_baits) is str : 
        bg_baits=[bg_baits] ; 

    #if len(cutoff) > 1 and len(cutoff) != len(bg_dses) : 
    #    sys.stderr.write('ERROR : If length of cutoff >=1, its length must match # of background dses.\n') ; 
    #    raise IndexError ; 

    #if len(query_dses) != len(query_baits) and len(query_dses) > 1 :
    #    sys.stderr.write('ERROR : Length of query_dses > 1 and query_baits must match.\n') ; 
    #    raise IndexError ; 
    #elif len(query_dses) != len(query_baits) : 
    #    # this is hackish in the extreme--
    #    # if you have ONE dataset (like a merged MS set) with which you want to interrogate multiple baits,
    #    # we make a list whose length matches that of the baits where each element points to the same data set
    #    tmp=query_dses[0] ; 
    #    query_dses=list() ; 
    #    while len(query_dses ) < len(query_baits) : 
    #        query_dses.append(tmp) ; 
    #if len(query_dses) != len(query_baits) :
    #    sys.stderr.write('ERROR : Length of bg_dses and bg_baits must match.\n') ; 
    #    raise IndexError ; 

    #sync lengths of query, query_baits
    if len(query_dses) != len(query_baits) : 
        if debug: 
            sys.stderr.write('Syncing query_dses (length {}) and query_baits (length {}).\n'\
             .format(len(query_dses),len(query_baits)))
        if len(query_dses) == 1 : 
            query_dses=[ query_dses[0] for i in mu.indices(query_baits) ] ; 
        elif len(query_baits) ==  1 : 
            query_baits=[ query_baits[0] for i in mu.indices(query_dses) ] ; 
        else : 
            sys.stderr.write('Can\'t sync lengths of query lists where both have len > 1.\n') ; 
            raise IndexError ; 

    #sync lengths of bg_dses, query_baits
    if len(bg_dses) != len(bg_baits) : 
        if debug: 
            sys.stderr.write('Syncing bg_dses (length {}) and bg_baits (length {}).\n'\
             .format(len(bg_dses),len(bg_baits)))
        if len(bg_dses) == 1 : 
            bg_dses=[ bg_dses[0] for i in mu.indices(bg_baits) ] ; 
        elif len(bg_baits) ==  1 : 
            bg_baits=[ bg_baits[0] for i in mu.indices(bg_dses) ] ; 
        else : 
            sys.stderr.write('Can\'t sync lengths of background lists where both have len > 1.\n') ; 
            raise IndexError ; 

    if len(query_dses) != len(cutoff) : 
        if debug: 
            sys.stderr.write('Syncing cutoffs (length {}) and controls (length {}).\n'\
             .format(len(bg_dses),len(bg_baits)))
        if len(cutoff ) == 1 : 
            cutoff=[ cutoff[0] for i in mu.indices(bg_baits) ] ; 
        else : 
            sys.stderr.write('Can\'t sync length of cutoffs to bg_dses where both have len > 1.\n') ; 
            raise IndexError ; 


    pseudos=[ np.mean([e.meanscore for e in bgd.nodes[bgb].edges['PSEUDO_00'] ]) for bgd,bgb in zip(bg_dses,bg_baits) ]  ;

    print([ query_dses[i].infilenames for i in range(0,len(query_dses)) ]) ; 

    if debug : 
        sew('PSEUDOCOUNTS :') ;  
        for ps in pseudos : 
            sew('\t{:6.2}'.format(ps)) ; 
        sew('\n') ; 
        

    out_pretup=list() ; 

    for qd,qb in zip(query_dses,query_baits) : 

        if as_dict : 
            outdict=dict() ; 
        else : 
            outset=set() ; 

        for ek,e in list(qd.edges.items()) : 

            if debug : 
                sew('EDGE : {} in {} against bait {}\n'.format(ek,query_dses.index(qd),qb)) ; 

            if qb not in {e.to.key,e.whence.key} : 
                if debug : 
                    sew('  N/A to current bait\n') ; 
                continue ; 
            elif e.to == e.whence  : 
                if selfs : 
                    nk = qb ; 
                else : 
                    continue ; 
            else :
                nk = next(iter( {e.to.key,e.whence.key} - {qb} )) ;
                if debug : 
                    sew('EDGE : {} in {} connects bait {} to {}\n'.format(ek,query_dses.index(qd),qb,nk)) ; 
                
                # whatever key is not the bait

            bg_sigs=list() ; 
            for bgd,bgb,bps in zip(bg_dses,bg_baits,pseudos) : 
                if nk in bgd.nodes[bgb].edges : 
                    bg_sigs.append(max([ bge.meanscore for bge in bgd.nodes[bgb].edges[nk] ])) ; 
                else:
                    bg_sigs.append(bps) ; 

            #!!
            #query_sig=qd.nodes[qb].edges[k].meanscore if qd.nodes[qb].edges.has_key(k) else 0.0 ; 
            query_sig=e.meanscore ;                            

            if debug : 
                sew('QUERY SIG : {:6.2} \n'.format(query_sig)) ; 
                sew('BG SIGS : ') ; 
                for s in bg_sigs : 
                    sew('\t{:6.2}'.format(s)) ; 
                sew('\n') ; 

            rd=[ ( query_sig - bgs ) / bgs for bgs in bg_sigs ]  ;
            if debug : 
                sew('REL DIFFS: ') ; 
                for s in rd : 
                    sew('\t{:6.2}'.format(s)) ; 
                sew('\n') ; 
            

            if as_dict : 
                outdict.update({ ek : rd }) ; 

            elif len(cutoff) == 1 :  
                # same cutoff for all bg sets
                for rdi in rd : 
                    if rdi < cutoff[0] : 
                        if debug : 
                            sew('    RD {:6.2} fails against cutoff {:6.2}\n'.format(rdi,cutoff[0])) ; 
                        break ; 
                else : 
                    if debug : 
                        sew('    Adding {}\n'.format(ek)) ; 
                    outset.add(ek) ;
            else : 
                # different cutoff for each bg set
                for rdi,cut in zip(rd,cutoff) : 
                    if rdi < cut : 
                        if debug : 
                            sew('    RD {:6.2} fails against cutoff {:6.2}\n'.format(rdi,cut)) ; 
                        break ; 
                else : 
                    if debug : 
                        sew('    Adding {}\n'.format(ek)) ; 
                    outset.add(ek) ; 

        if as_dict : 
            out_pretup.append( outdict ) ;
        else : 
            #print type(outset) ;
            out_pretup.append( outset ) ;

    return tuple(out_pretup) ; 


def zfilter(query_ds,bg_dses,query_bait,bg_bait,logp = -3,as_dict=False,percentile=50,\
             directed=None,qual=None,debug=False,logzero=1e-300) : 
    """
        This works differently than bgsubtract.
        You get ONE query ds and ONE query bait and ONE cutoff.
        Your bg_dses are pooled and the z xcore of your query signals are measured against
        ONE bait.

        As usual, having as_dict as False returns a list of keys that cleared the filter ; 
        having it as True ignores the Z filter and outputs a dict as such : 
            key : z
    """

    from scipy.stats.distributions import truncnorm

    dbgstr='DEBUG> interactors_extras.zfilter :'

    # because we have to do this by edge key rather than node key
    qeks=[ e.key for e in list(query_ds.edges.values())\
             if query_bait in {e.to.key,e.whence.key} and\
             ( directed is None or ( e.directed == directed and e.whence.key == query_bait )) and\
             ( qual     == None   or e.qual     == qual ) ] ; 

    if len(qeks) == 0  :
        sys.stderr.write('interactors_extras.zfilter : Length of fetched edge values is 0!\n'\
                         'Dumping states of relevant variables.\n') ; 
        sys.stderr.write('    query_ds  :   '+repr(query_ds)+'\n') ; 
        sys.stderr.write('    bg_dses   :   '+repr(bg_dses)+'\n') ; 
        sys.stderr.write('    query_bait:   '+repr(query_bait)+'\n') ; 
        sys.stderr.write('    bg_bait   :   '+repr(bg_bait)+'\n') ; 
        sys.stderr.write('    directed  :   '+repr(directed)+'\n') ; 
        sys.stderr.write('    qual      :   '+repr(qual   )+'\n') ; 
        return

    signals=dict() ; 
    preys=dict() ;
    ctrl_means=dict() ; 
    ctrl_variances=dict() ; 

    for ek in qeks : 
        signals.update({ ek : query_ds.edges[ek].meanscore }) ; 
        preys.update({ ek : query_ds.edges[ek].connects_key(query_bait) }) ; 

    for ek in qeks : 

        prey=preys[ek] ; 

        ctrl_vals=[ sum([ e.meanscore for e in d.nodes[bg_bait].edges[prey] ]) if prey in d.nodes else \
                    0.0 for d in bg_dses] ; 
        # pseudocounts were found to fit more poorly

        ctrl_means.update({     ek : np.mean(ctrl_vals) }) ; 
        ctrl_variances.update({ ek :  np.var(ctrl_vals) if any([ x > 0.0 for x in ctrl_vals ]) else None }) ; 

    pseudomean  = np.mean([ e.meanscore for d in bg_dses for e in d.nodes[bg_bait].edges['PSEUDO_00'] ]) ; 
    pseudovar   = np.var([  e.meanscore for d in bg_dses for e in d.nodes[bg_bait].edges['PSEUDO_00'] ]) ;

    min_var     =min([ v for v in ctrl_variances.values() if v is not None ]) ; 

    if debug : 
        sys.stderr.write(dbgstr+'pseudomean {} pseudovar {} min_var {}\n'.\
            format(pseudomean,pseudovar,min_var)) ; 

    # pass to account for all-pseudocount observations
    for ek in set(qeks) - {'PSEUDO_00'} : 
        if ctrl_means[ek] == 0.0 and ctrl_variances[ek] == None : 
            # comment later
            #print 'Edge ',ek,'connects to prey unobserved in controls.'
            ctrl_variances.update({ ek : pseudovar }) ; 

    Z   = lambda k : ( signals[k] - ctrl_means[k] ) / np.sqrt(ctrl_variances[k]) ; 
    Zto0= lambda k : -1.0 * (ctrl_means[k]) / np.sqrt(ctrl_variances[k]) ; 

    if as_dict : 

        out=dict() ;
        for ek in set(qeks) - {'PSEUDO_00'} : 
           out.update({ ek : np.log10(max([1.0 - truncnorm.cdf(Z(ek),Zto0(ek),np.inf),logzero])) }) ; 
        return out ; 

    else : 
        out=set() ; 

        for ek in set(qeks) - {'PSEUDO_00'} : 
            if np.log10(max([1.0 - truncnorm.cdf(Z(ek),Zto0(ek),np.inf),logzero])) <= logp : 
                out.add(ek) ; 

    refuse=set(qeks) - out ; 
    if percentile > 0  and len(refuse) > 0 : 
        wheat=set() ; # to be separated from chaff
        #thresh = np.median([ signals[ek] for ek in refuse ]) ; 
        thresh = np.percentile([ signals[k] for k in refuse ],percentile) ; 
        if debug : 
            sys.stderr.write(dbgstr+'l(median) of rejects is {}\n'.format(thresh)) ; 
        #print thresh

        wheat={ ek for ek in out if signals[ek] > thresh }

        return wheat ; 

    else : 
        return out

def print_springs( edges, fname = "", print_headers = True, sep = '\t', print_weights = False,
                   print_total_scores = False, print_mean_scores = False, print_organisms = False,
                   print_source = False, inter_string = 'pp', transform_scores = None,
                   print_quals = True, print_pps = False, print_saint = False, print_bkg = False ):

    def spring_transform(edge) : 

        if edge.qual.lower() in {'em', 'emili', 'complexes'} : 
            #print( edge.qual.lower() + '10')
            return 10  ;
        elif edge.qual.lower() in {'biogrid', 'bg', 'bp', 'bioplex', 'hippie'}:
            #print( edge.qual.lower() + '1')            
            return 1
        else:
            #print(str(8.5+round(np.log10(edge.meanscore), 2)))
            return 8.5 + round(np.log10(edge.meanscore), 2) 

    if ( not fname ) :
        f = sys.stdout
    else :
        f = open(fname,"w") ; 

    if ( print_headers) :
        f.write("Left{}Inter{}Right".format(sep,sep)) ;
        if ( print_quals ) :
            f.write("{}Qual".format(sep)) ;
        if ( print_weights ) :
            f.write("{}Wt".format(sep)) ;

        # guaranteed mean printing
        f.write("{}Mean".format(sep)) ;

        if ( print_total_scores) :
            f.write("{}Total".format(sep)) ;
        if ( print_organisms ) :
            f.write("{}OrgA{}OrgB".format(sep,sep)) ; 
        if ( print_source ) :
            f.write("{}Source".format(sep,sep)) ; 
        if print_pps :
            if print_saint:
                f.write("{}fdr".format(sep))
            f.write("{}logp".format(sep)) ; 

        # guaranteed springiness printing
        f.write("{}Spring".format(sep)) ;

        if print_saint:
            f.write("{}avgP".format(sep)) ;

        if print_bkg:
            f.write("{}Background Values".format(sep))
            f.write("{}pVal".format(sep))            
            f.write("{}Spectral Counts".format(sep))
            f.write("{}Coverage".format(sep))
            f.write("{}Unique peptides".format(sep))
            
            
        f.write("\n") ;

    for edge in edges : 

        #edget=tuple(edge.nodes) ;
        thedir = '>' if edge.directed else '^'

        f.write("{}{}{}{}{}{}".format(edge.whence.official,sep,thedir,inter_string,sep,edge.to.official)) ;

        if ( print_quals ):
            f.write("{}{}".format(sep,edge.qual)) ;
        if ( print_weights ):
            f.write("{}{}".format(sep,edge.weight)) ;

        # guaranteed mean printing
        if transform_scores : 
            f.write("{}{}".format(sep,transform_scores(edge.meanscore))) ;
        else : 
            f.write("{}{}".format(sep,edge.meanscore)) ;

        if ( print_total_scores):
            if transform_scores : 
                f.write("{}{}".format(sep,transform_scores(edge.totalscore))) ;
            else : 
                f.write("{}{}".format(sep,edge.totalscore)) ;

        if ( print_organisms):
            if ( len(edge.nodes) == 2 ):
                f.write("{}{}{}{}".format(sep,edget[0].organism,sep,edget[1].organism)) ;
            else :
                f.write("{}{}{}{}".format(sep,edget[0].organism,sep,edget[0].organism)) ;
        if ( print_source):
            f.write("{}{}".format(sep,edge.source)) ;
        if print_pps :
            if print_saint:
                f.write("{}{}".format(sep,edge.p))
            f.write("{}{}".format(sep,-1*np.log10(edge.p)))

        # guaranteed springiness printing

        f.write("{}{}".format(sep,spring_transform(edge))) ; 

        if print_saint:
            f.write("{}{}".format(sep,edge.avgP))

        if print_bkg:
            f.write("{}{}".format(sep,edge.bkg))
            f.write("{}{}".format(sep,edge.p_ori))
            edge.spc   =  '|'.join([ re.sub( r'^.+raw_(\d+)_?.*$', r'\1', i.tags) for i in edge.interactions ])
            edge.cov   =  '|'.join([ re.sub( r'^.+cov_([\d\.]+)_?.*$', r'\1', i.tags) for i in edge.interactions ])
            edge.upept =  '|'.join([ re.sub( r'^.+upept_(\d+)_?.*$', r'\1', i.tags) for i in edge.interactions ])            
            f.write("{}{}".format(sep,edge.spc))
            f.write("{}{}".format(sep,edge.cov))
            f.write("{}{}".format(sep,edge.upept))            
            
        f.write("\n") ;

    f.close()

    
def dombuddies(infile,keys,debug=False,rooted=False) : 
    warn(DeprecationWarning('I (Mark) can\'t vouch for dombuddies right now 3/1/2016') )
    #import rbase
    from lib.markutils import b4us,afterus
    rbase.load('cdd')
    rbase.load('hmg') ; 

    myds=I.dataSet() ; 

    buddies=list() ; 
    alldoms=set() ;
    kroots=dict() ; 
    # maps domains to keys whose genes contain that domain
    
    for line in infile : 
        linel=line.strip().split('\t') ;

        buddies.append(( linel[0],linel[1] )) ; 

        alldoms.update({ linel[0],linel[1] }) ; 
        kroots.update({ linel[0] : set() }) ; 
        kroots.update({ linel[1] : set() }) ; 

    for k in keys : 

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # add node for this key to dataset
        sym=b4us(k) ; 
        eid=afterus(k) ; 

        org=rbase.hmg['EID'][eid]['Taxon'] if eid in rbase.hmg['EID'] else '00' ;

        myds.nodes.update({ k : I.node(entrez=eid,official=sym,\
                                 organism=org,key=k) }) ; 
        if debug : 
            sys.stdderr.write('DEBUG: Added node with key '+k+'\n') ; 

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # add dict of roots for this dataest
        #kroots.update({ k : { rbase.cdd['ByPssm'][c]['Root']['Pssm']\
        #                      for c in rbase.hmg['EID'][eid]['CDD'] if\
        #                      c in rbase.cdd['ByPssm'] } }) ; 
        if eid not in rbase.hmg['EID'] or rbase.hmg['EID'][eid]['CDD'] is None : 
            continue ; 


        for c in rbase.hmg['EID'][eid]['CDD'] : 
            if rooted : 

                if c in rbase.cdd['ByPssm'] and rbase.cdd['ByPssm'][c]['Root']['Pssm'] in kroots : 
                    theroot=rbase.cdd['ByPssm'][c]['Root']['Pssm'] ; 

                    kroots[theroot].add(k) ; 

            else : 
                if c in kroots : 
                    kroots[c].add(k) ; 
                

    for d1,d2 in buddies : 

        j=0 ; 
        for k1 in kroots[d1] : 
            for k2 in kroots[d2] : 
                
                thisinteraction=I.interaction(myds.nodes[k1],myds.nodes[k2], \
                 interID    =   'dp_'+d1+'^'+d2+'_'+repr(j),\
                 qualifications = 'dp') ; 

                myds.the_data.update({ thisinteraction.interID : thisinteraction }) ; 

                j += 1 ;

                if thisinteraction.edgekey() not in myds.edges :
                # no need to check for inversion

                    newedge=I.bgedge(interaction=thisinteraction,\
                            directed=thisinteraction.directed,\
                            qual=thisinteraction.qualifications) 

                    myds.edges.update({ newedge.key : newedge }) ; 

                    newedge.whence.partners.update({ newedge.to.key : newedge.to }) ; 
                    newedge.to.partners.update({ newedge.whence.key : newedge.whence }) ; 

                    if newedge.whence.key not in newedge.to.edges : 
                        newedge.to.edges.update({ newedge.whence.key : { newedge } }) ; 
                    else : 
                        newedge.to.edges[newedge.whence.key].add(newedge) ; 

                    if newedge.to.key not in newedge.whence.edges : 
                        newedge.whence.edges.update({ newedge.to.key : { newedge } }) ; 
                    else : 
                        newedge.whence.edges[newedge.to.key].add(newedge) ; 
                
                else :
                    myds.edges[thisinteraction.edgekey()].add_interaction(thisinteraction) ; 


                thisinteraction.nodeA.interactions.update({ thisinteraction.interID : thisinteraction }) ; 
                thisinteraction.nodeB.interactions.update({ thisinteraction.interID : thisinteraction }) ;

    infile.seek(0) ; 

    return myds

   #def cliquebait(nodes,wes=None,mindegree=2) :

   #    cliques=list()
   #    # TODO : implement bron-kerbosch with a given set of edges

   #    def bronkerbosch(output,totry,excluded) : 

   #        if not totry and not excluded : 
   #            return output : 
   #        else : 
   #            for v in set(totry) : 

   #                # TODO implement function that sees if a node is a partner
   #                # within a given edge set
   #                # and while you're at it, one that gives the score from one node to another
   #                # within a given edge set.
   #                # the function we ACTUALLY need should return v's neighbors along a given set of edges

   #                output = bronkerbosch( output | {v}, totry & vn , excluded & vn ) ; 

   #                totry  -= v  ;
   #                excluded += v ; 

def zfilter_mp(query_ds,bg_dses,query_bait,bg_bait,logp = -3,as_dict=False,percentile=50,\
             directed=None,qual=None,debug=False,logzero=1e-300) : 
    """
        This works differently than bgsubtract.
        You get ONE query ds and ONE query bait and ONE cutoff.
        Your bg_dses are pooled and the z xcore of your query signals are measured against
        ONE bait.

        As usual, having as_dict as False returns a list of keys that cleared the filter ; 
        having it as True ignores the Z filter and outputs a dict as such : 
            key : z
    """
    warn(DeprecationWarning('I (Mark) can\'t vouch for zfilter_mp right now 3/1/2016') )

    #mplog10 =  lambda x : mpmns.log(x)/mpmns.log(10)
    #mptndpf  = lambda z : mpmns.npdf(x) * ( 1 + 1/mpmns.ncdf(
    def mp_tnpdf(z,z0) : 
        if z < z0 : 
            return mp.ninf ; 

        else : 
            return mpmns.npdf(z) / ( mpmns.ncdf(mpmns.inf) - mpmns.ncdf(z0) ) ;

    def mp_tncdf(z,z0) :
        if z < z0 : 
            return mp.ninf ; 
        else : 
            return ( mpmns.ncdf(z) - mpmns.ncdf(z0) )  / ( mpmns.ncdf(mpmns.inf) - mpmns.ncdf(z0) ) ;

    dbgstr='DEBUG> interactors_extras.zfilter :'

    # because we have to do this by edge key rather than node key
    qeks=[ e.key for e in list(query_ds.edges.values())\
             if query_bait in {e.to.key,e.whence.key} and\
             ( directed is None or ( e.directed == directed and e.whence.key == query_bait )) and\
             ( qual     == None   or e.qual     == qual ) ] ; 

       #if len(qeks) == 0  :
       #    print query_ds
       #    print bg_dses
       #    print query_bait
       #    print bg_bait
       #    print z
       #    print directed
       #    print qual
       #    return

    signals=dict() ; 
    preys=dict() ;
    ctrl_means=dict() ; 
    ctrl_variances=dict() ; 

    for ek in qeks : 
        signals.update({ ek : query_ds.edges[ek].meanscore }) ; 
        preys.update({ ek : query_ds.edges[ek].connects_key(query_bait) }) ; 

    for ek in qeks : 

        prey=preys[ek] ; 

        ctrl_vals=[ mpmns.mpf(sum([ e.meanscore for e in d.nodes[bg_bait].edges[prey] ])) if prey in d.nodes else \
                    0.0 for d in bg_dses] ; 
        # pseudocounts were found to fit more poorly

        ctrl_means.update({     ek : mpmns.mpf(np.mean(ctrl_vals)) }) ; 
        ctrl_variances.update({ ek :  mpmns.mpf(np.var(ctrl_vals)) if any([ x > 0.0 for x in ctrl_vals ]) else None }) ; 

    pseudomean  = np.mean([ mpmns.mpf(e.meanscore) for d in bg_dses for e in d.nodes[bg_bait].edges['PSEUDO_00'] ]) ; 
    pseudovar   = np.var([  mpmns.mpf(e.meanscore) for d in bg_dses for e in d.nodes[bg_bait].edges['PSEUDO_00'] ]) ;

    min_var     =min(ctrl_variances.values()) ; 

    if debug : 
        sys.stderr.write(dbgstr+'pseudomean {} pseudovar {} min_var {}\n'.\
            format(pseudomean,pseudovar,min_var)) ; 

    # pass to account for all-pseudocount observations
    for ek in set(qeks) - {'PSEUDO_00'} : 
        if ctrl_means[ek] == 0.0 and ctrl_variances[ek] == None : 
            # comment later
            #print 'Edge ',ek,'connects to prey unobserved in controls.'
            ctrl_variances.update({ ek : pseudovar }) ; 

    Z   = lambda k : ( signals[k] - ctrl_means[k] ) / np.sqrt(ctrl_variances[k]) ; 
    Zto0= lambda k : -1.0 * (ctrl_means[k]) / np.sqrt(ctrl_variances[k]) ; 

    if as_dict : 

        out=dict() ;
        for ek in set(qeks) - {'PSEUDO_00'} : 
           out.update({ ek : mpmns.log10(1.0 - mp_tncdf(Z(ek),Zto0(ek))) }) ; 
        return out ; 

    else : 
        out=set() ; 

        for ek in set(qeks) - {'PSEUDO_00'} : 
            if mpmns.log10(1.0 - mp_tncdf(Z(ek),Zto0(ek))) <= logp : 
                out.add(ek) ; 

    if percentile > 0 : 
        refuse=set(qeks) - out ; 
        wheat=set() ; # to be separated from chaff
        #thresh = np.median([ signals[ek] for ek in refuse ]) ; 
        thresh = np.percentile([ signals[k] for k in refuse ],percentile) ; 
        if debug : 
            sys.stderr.write(dbgstr+'l(median) of rejects is {}\n'.format(thresh)) ; 
        #print thresh

        wheat={ ek for ek in out if signals[ek] > thresh }

        return wheat ; 

    else : 
        return out

def madfilter(dataset,ctrl_fname,baitkey,qual=None,directed=False,as_dict=False,alpha=0.05,debug=False) : 
    import pickle

    import os
    from scipy.stats.distributions import norm
    if not os.path.isfile(ctrl_fname) and not \
     os.path.isfile( cf.controlFiles + ctrl_fname) :
        raise FileNotFoundError('Could not find '+ctrl_fname) ;
    elif os.path.isfile( cf.controlFiles + ctrl_fname) : 
        fn = cf.controlFiles + ctrl_fname ; 
    else : 
        fn=ctrl_fname ;

    if debug : 
        sys.stderr.write('DEBUG> interactors_extras.madfilter : fname is '+fn+'\n') ;


    with open(fn,'rb') as f : 
        (syms,medians,mads)=pickle.load(f) ;

    if debug :
        sys.stderr.write('DEBUG> interactors_extras.madfilter :\n'+\
         '        syms: '+str(len(syms))+' medians: '+str(len(medians))+' mads: '+str(len(mads))+\
         '\n') ;

    symset=set(syms) ; 

    pseudoindex=syms.index('PSEUDO') ;
    maddict=dict() ; 
    peedict=dict() ; 
    if not as_dict : 
        outedges=set() ;

    for e in { e for es in dataset.nodes[baitkey].edges.values() for e in es  } : 
        if qual and e.qual != qual : 
            if debug : 
                sys.stderr.write('DEBUG> interactors_extras.madfilter : edge '+e.key+\
                ' skipped due to invalid qualifier\n'+e.qual) ;
            continue ; 
        if directed and e.to.key == baitkey : 
            if debug : 
                sys.stderr.write('DEBUG> interactors_extras.madfilter : edge '+e.key+\
                ' skipped due to wrong direction\n') ;
            continue ;

        if e.to.official not in symset : 
            madscore= ( np.log10(e.meanscore) - medians[pseudoindex] )/ mads[pseudoindex] ; 
            if debug : 
                sys.stderr.write('DEBUG> interactors_extras.madfilter : edge '+e.key+\
                'given score '+'{:8.6}'.format(madscore)+' (pseudocounted) \n') ;
        elif e.to.official in symset : 
            i=syms.index(e.to.official)
            madscore= ( np.log10(e.meanscore) - medians[i] )  / mads[i] ; 
            if debug : 
                sys.stderr.write('DEBUG> interactors_extras.madfilter : edge '+e.key+\
                'given score '+'{:8.6}'.format(madscore)+'  \n') ;

        maddict.update({ e.key : madscore }) ;
        peedict.update({ e.key : 1-norm.cdf(madscore) })
        #elif madscore > sds * 1.48 : 
            #outedges.add( e.key ) 

    eks=list(maddict.keys()) ; 
    mads=list(maddict.values()) ; 
    pees=list(peedict.values()) ; 
    from statsmodels.sandbox.stats.multicomp import multipletests

    rejects,adjpees=multipletests(pees,alpha=alpha)[0:2] ;

    if as_dict : 
        apdict=dict(zip(eks,adjpees)) ; 
        return (maddict,apdict) ;
    else :
        return { eks[x] for x in range(len(eks)) if rejects[x] } ; 

def dataset_edges_for_bait( ds, baitkey, qual, directed, debug = False ):

    ek_ps = dict() # edgekey --> preysymbol
    ek_ms = dict() # edgekey --> meanscore
    ps_ms = dict() # preysymbol --> meanscore

    for e in { e for es in ds.nodes[baitkey].edges.values() for e in es  } :
        if qual and e.qual != qual : 
            if debug : 
                sys.stderr.write('DEBUG> interactors_extras.madfilter_corr : edge '+e.key+\
                                 ' skipped due to invalid qualifier: '+qual + '\n') ;
            continue ; 
        elif directed and e.to.key == baitkey : 
            if debug : 
                sys.stderr.write('DEBUG> interactors_extras.madfilter_corr : edge '+e.key+\
                ' skipped due to wrong direction\n') ;
            continue ;
        else:
            # for the reason behind this constant, cf
            # https://en.wikipedia.org/wiki/Median_absolute_deviation#Relation_to_standard_deviation
            if debug : 
                sys.stderr.write( 'keep this edge: ' + e.key + '\n')
            ek_ps.update({ e.key : e.to.official })
            ek_ms.update({ e.key : np.log10(e.meanscore) }) ;
            ps_ms.update({ e.to.official : np.log10(e.meanscore) })

    return ( ek_ps, ek_ms, ps_ms )

def mad(series) : 
    return np.percentile(np.abs(series-np.percentile(series,50)),50) ;

def get_ortholog( symbol, conv, targ ):

    conv = { conv[k][0]: conv[k][1] for k in conv if conv[k][0] is not None }
    
    if symbol in conv:
        if len(conv[symbol]) == 1:
            # if there is a single match, use that
            symbol = conv[symbol]
        else:
            for s in conv[symbol]:
                # select the one with the same name
                if s.lower() == symbol.lower():
                    symbol = s
                    break
            else:
                # else, select the one with the smallest eid
                eid = str(sorted([int(targ['Symbol'][k]['EID']) for k in conv[symbol]])[0])
                symbol = targ['EID'][eid]['Symbol']
                
    return symbol
    

def read_control_data( ctrl_fname, convert, debug  ):
    # read in control file
    
    if not os.path.isfile(ctrl_fname) and not os.path.isfile(CONTROL_FILES + ctrl_fname) :
        raise FileNotFoundError('Could not find '+ctrl_fname) ;
    elif os.path.isfile(CONTROL_FILES + ctrl_fname) : 
        fn = CONTROL_FILES + ctrl_fname ; 
    else : 
        fn = ctrl_fname ;

    # if the reference file is just a list of .i files, create the pickled datafile
    if re.match( '.*\.txt$', fn ):
        fn = createReferenceFile( fn )

    with open(fn,'rb') as f : 
        ( allsyms, allfns, loggrid ) = pickle.load(f) ;

    if debug :
        sys.stderr.write( 'DEBUG> interactors_extras.madfilter_corr : control fname is ' + fn + '\n' )
        sys.stderr.write('DEBUG> interactors_extras.madfilter_corr : control\n'+
                         '        syms: '+str(len(allsyms))+' files: '+str(len(allfns))+'\n') ;

    if convert is not None:
        asfile = fn + '_allsymb'
        if os.path.isfile( asfile ):
            with open(asfile, 'rb') as x:
                allsymb = pickle.load(x)
        else:
            rbase.load('m2hs') if convert == 'm2h' else rbase.load('h2ms')
            conv   = rbase.m2hs if convert == 'm2h' else rbase.h2ms
            rbase.load('hsg')  if convert == 'm2h' else rbase.load('mmg')
            targ   = rbase.hsg  if convert == 'm2h' else rbase.mmg
            allsym = [get_ortholog( symb, conv, targ ) for symb in allsyms ]
            with open(asfile, 'wb') as x:
                pickle.dump(allsym, x)
        
    return( allsyms, allfns, loggrid )
        
    
def madfilter_corr( dataset,                # network dataset to process, interactors.dataSet instance
                    ctrl_fname,             # control file name to use, pickled syms,medians mads (*.cp2))
                    #                         or a *.txt file listing one ifile per line 
                    baitkey,                # key for the bait in dataset
                    qual     = None,        # qualifier for this bait,
                    directed = False,       # is the network directed/non-
                    as_dict  = False,       # returns two dictionary objects of mad and p-value for each edge
                    alpha    = 0.05,        # fwer or fdr deepending on method
                    convert  = None,        # whether to convert the control file to another org 
                    debug    = False,       #
                    maxcorr  = 0.75,        #
                    method   = 'fdr_bh',    # method for multiple hypothesis testing correction
                    assign_edge_ps  = True  # whether or not to modify 'p' attribute of tested edges
) : 

    """
        madfilter_corr(dataset,ctrl_fname,baitkey,qual=None,directed=False,\
            as_dict=False,alpha=0.05,debug=False,maxcorr=0.75)
        dataset : interactors.dataSet instance
        ctrl_fname : pickled syms,medians mads 
        baitkey : key of bait
    """

    ( allsyms, allfns, loggrid ) = read_control_data( ctrl_fname, convert, debug )
    symset      = set(allsyms) ; 
    pseudoindex = allsyms.index('PSEUDO') ;
    #pseudoscore = dataset.nodes[baitkey].edges.get('PSEUDO_00', PSEUDO_DEFAULT)
    idstr       = baitkey + '>' + 'PSEUDO_00' + ':' + qual
    pseudoscore = dataset.edges[idstr].meanscore if idstr in dataset.edges else PSEUDO_DEFAULT
    pseudoscore = np.log10( pseudoscore )
    
    if not as_dict : 
        outedges = set() ;

    # this loop now does X things : 
    #   1) get a vector of scores to test
    #       -specified by bait, direction, qual
    #   2) put those scores into dicts :

    ek_ps, ek_ms, ps_ms = dataset_edges_for_bait( dataset, baitkey, qual, directed, debug )

    # NEXT : remove datasets from background that are highly correlated with the query dataset
    v = np.array([ ps_ms.get( allsyms[x], ps_ms['PSEUDO']) for x in range(len(allsyms)) ])

    highly_correlated_fn_indices = [ x for x in range(len(allfns)) if np.corrcoef(v,loggrid[:,x])[0,1] > maxcorr ]

    if debug:
        corrs = [str(np.corrcoef(v,loggrid[:,x])[0,1]) for x in range(len(allfns))]
        sys.stderr.write('interactors_extras.madfilter_corr : the following correlations '
                         + 'were observed: \n' )
        for i in range(len(allfns)):
                         sys.stderr.write( '        ' + allfns[i] + "\t" + corrs[i] + '\n' )

    if highly_correlated_fn_indices : 
        sys.stderr.write('interactors_extras.madfilter_corr : the following experiments\n'+\
        '        are too highly correlated with the current one for inclusion\n') ;
        for x in highly_correlated_fn_indices  : 
            sys.stderr.write('          '+allfns[x]+'\n') ;

    # remove the highly correlated expts
    compgrid = np.delete(loggrid,highly_correlated_fn_indices,1)

    # calculate mad scores and p values
    maddict  = dict()
    peedict  = dict()
    bkgdict  = dict()

    for ek in ek_ps.keys() :
        e = dataset.edges[ek] ;
        i = pseudoindex
        #if e.to.official not in symset : 
        #    madscore = (ek_ms[ek]- np.median(compgrid[pseudoindex,:])) / mad(compgrid[pseudoindex,:])/ 1.48 ;
        #    bkgvls   = compgrid[pseudoindex,:]
        #    if debug :
        #        sys.stderr.write('DEBUG> interactors_extras.madfilter_corr : edge '+e.key+\
        #        ' not in background - given score '+'{:8.6}'.format(madscore)+' (pseudocounted) \n') ;

        #elif e.to.official in symset :
        if e.to.official in symset : 
            i    = allsyms.index( e.to.official )

        madscore = ( ek_ms[ek] - np.median(compgrid[i,:])) / mad(compgrid[i,:]) / 1.48 ;
        p05_cut  = 1.65 * 1.48 * mad(compgrid[i,:]) + np.median(compgrid[i,:])
        bkgvls   = compgrid[i,:]
        pval     = 1-norm.cdf(madscore)

        maddict.update({ e.key : madscore }) ;
        peedict.update({ e.key : pval })
        bkgdict.update({ e.key : "|".join( map( str, bkgvls)) })

        if debug : 
            sys.stderr.write('DEBUG> interactors_extras.madfilter_corr : edge='+e.key +
                             ' madscore='+'{:8.6}'.format(madscore)+'; nsaf=' +
                             '{:8.6}'.format(ek_ms[ek]) + ' p05_cut='+ '{:8.6}'.format(p05_cut) + ' pval=' +
                             '{:8.6}'.format(pval) + '  \n')
            
        
    eks  = list( maddict.keys() ) 
    mads = list( maddict.values() ) 
    pees = list( peedict.values() )
    bkgs = list( bkgdict.values() )

    rejects, adjpees = multipletests( pees, alpha = alpha, method = method )[0:2]

    passed = { eks[x] for x in range(len(eks)) if rejects[x] }

    if assign_edge_ps : 
        for i in range(len(eks)) : 
            dataset.edges[eks[i]].p     = adjpees[i]
            dataset.edges[eks[i]].p_ori = pees[i]
            dataset.edges[eks[i]].bkg   = bkgs[i]
            
    #adjval = { pees[x] for x in range(len(eks)) if rejects[x] }
    if debug : 
        import colorama
        sys.stderr.write('DEBUG> Results summary:\n') ;
        print('index','edge_key','p_value','std','rejected','p_adj','pNSAF',sep='\t')
        for x in range(len(eks)) : 
            e = dataset.edges[eks[x]]
            if rejects[x] : 
                print(colorama.Back.BLUE+repr(x),eks[x],peedict[eks[x]],maddict[eks[x]],rejects[x],\
                 adjpees[x],repr(-1*np.log10(e.meanscore))+colorama.Back.BLACK,sep='\t')
            else : 
                #print(colorama.Back.BLUE,end='')
                print(x,eks[x],peedict[eks[x]],maddict[eks[x]],rejects[x],\
                 adjpees[x],-1*np.log10(e.meanscore),sep='\t')
                #print(colorama.Style.RESET_ALL,end='')

    if as_dict : 
        apdict = dict(zip(eks,adjpees)) ; 
        return (maddict,apdict) ;
    else :
        return passed

def createReferenceFile( inputFile, overWrite = False ):

    if re.match( '.*cp2$', inputFile ):
        outf        = inputFile
    else :
        outf        = re.sub( 'txt$', 'cp2', inputFile )

    if not isfile( outf ) or overWrite :
        
        alldatadict = dict()
    
        with open( inputFile ) as f : 
            for line in f : 
                if re.match( '^#.*|^\s*$', line ):
                    continue
                alldatadict.update( {line.strip():dict()} ) ;

        for fn in alldatadict.keys() :
            f = open( cf.ifilesPath + fn ) 
            f.readline()
            for ll in mu.tabulate(f) : 
                alldatadict[fn].update({ ll[2] : np.float(ll[3]) })
            f.close() ;

        pseudodict     = { k : alldatadict[k]['PSEUDO'] for k in alldatadict }
        pskeys         = list(pseudodict.keys())
        pslogvals      = np.log10(list(pseudodict.values()))
        pslogmad       = mad(pslogvals) ; 
        pslogmedian    = np.percentile(pslogvals,50)
        pslvps_hi      = 1-norm.cdf((pslogvals-pslogmedian)/pslogmad)
        #pslvps_lo=1-norm.cdf((pslogvals-pslogmedian)/pslogmad)

        rejected_ds_hi = multipletests(pslvps_hi,alpha=0.05)[0]
        #rejected_ds_lo=multipletests(pslvps_lo,alpha=0.05)[0]

        allfns         = list()
        allsyms        = set()
        for i in range(len(rejected_ds_hi)) : 
            if rejected_ds_hi[i] : 
                sys.stderr.write('WARNING: '+pskeys[i]+' is a pseudocount outlier '+\
                                 'and as such will NOT be incorporated\n')
            else : 
                allfns.append(pskeys[i]) ; 
                allsyms |= set(alldatadict[pskeys[i]].keys())

        allsyms        = list(allsyms) ; 

        # creation of the all-symbols all-fns grid
        #allsyms=list({ k for fn in alldatadict.keys() for k in alldatadict[fn].keys() })
        #allfns=list(alldatadict.keys())

        grid           = np.zeros((len(allsyms),len(allfns)))

        for i in range(len(allsyms)) : 
            for j in range(len(allfns)) : 
                grid[i][j] = alldatadict[allfns[j]].get(allsyms[i],alldatadict[allfns[j]]['PSEUDO'])

        loggrid        = np.log10(grid) ;
        logmads        = np.fromiter( ( mad(loggrid[i,:]) for i in range(len(allsyms)) ),dtype=np.float) ; 
        logmeans       = np.fromiter( ( np.mean(loggrid[i,:]) for i in range(len(allsyms)) ),dtype=np.float) ; 
        logmedians     = np.fromiter( ( np.percentile(loggrid[i,:],50) for i in range(len(allsyms)) ),dtype=np.float) ; 

        with open( outf,'wb') as f : 
            pickle.dump( (allsyms, allfns, loggrid), f ) ; 

    return( outf )

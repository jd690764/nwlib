import pickle
import readline
import colorama
import sys
import os

def arm(loc) : 
    out=''

    if loc is None : 
        return None ; 

    for c in loc : 
        try :
            int(c) ;    
            out += c; 
        except : 
            out +=c; 
            return out ; 
    else : 
        return out ; 

def chr(loc) : 
    out=''

    if loc is None : 
        return None ; 

    for c in loc : 
        try :
            int(c) ;    
            out += c; 
        except : 
            return out ; 
    else : 
        return out ; 

refsuite=None ; 
biogrid=None ; 
preppi=None ; 
dompairs=None ; 

OUTCOLOR=colorama.Fore.RED
OUTSTYLE=colorama.Style.BRIGHT

#REFERENCEPATH='/mnt/reference/' ;
REFERENCEPATH='/usr/local/share/py/djscripts/data/pickles/' ;

FILESDICT={\
'hsg' : REFERENCEPATH + 'hsg_latest' ,\
'mmg' : REFERENCEPATH + 'mmg_latest' ,\
'hmg' : REFERENCEPATH + 'hsmmg_latest' ,\
'hsp' : REFERENCEPATH + 'hsp_latest' ,\
'mmp' : REFERENCEPATH + 'mmp_latest' ,\
'cdd' : REFERENCEPATH + 'cdd_latest' ,\
'h2m' : REFERENCEPATH + 'h2m_latest' ,\
'm2h' : REFERENCEPATH + 'm2h_latest' ,\
'dup' : REFERENCEPATH + 'dup_latest' ,\
}

hsg=None ; 
mmg=None ; 
hmg=None ; 
hsp=None ; 
mmp=None ; 
cdd=None ; 
m2h=None ; 
h2m=None ; 
dup=None ; 
others=dict() ;  

ATTRDICT={\
'hsg' : hsg,\
'mmg' : mmg,\
'hmg' : hmg,\
'hsp' : hsp,\
'mmp' : mmp,\
'cdd' : cdd,\
'm2h' : m2h,\
'h2m' : h2m,\
'dup' : dup,\
'others' : others,\
}

VALID_CDDS=set() ; 

def refroots() : 

    global hsg ; 
    global mmg ; 
    global hmg ;
    global hsp ; 
    global mmp ; 
    global cdd ; 
    global m2h ; 
    global h2m ;
    global dup ; 


    hsg=ATTRDICT['hsg'] ; 
    mmg=ATTRDICT['mmg'] ; 
    hmg=ATTRDICT['hmg'] ; 
    hsp=ATTRDICT['hsp'] ; 
    mmp=ATTRDICT['mmp'] ; 
    cdd=ATTRDICT['cdd'] ; 
    m2h=ATTRDICT['m2h'] ; 
    h2m=ATTRDICT['h2m'] ; 
    dup=ATTRDICT['dup'] ; 
    others=ATTRDICT['others'] ; 

        

def load(db) : 

    from os.path import isfile

    try : 
        if db not in ATTRDICT and isfile(db) : 
            fobj=open(db,'rb') ; 
            ATTRDICT['others'].update({ fobj.name.split('/')[-1] : pickle.load(fobj) }) ; 
            fobj.close() ; 
        elif ATTRDICT[db] is None : 
            fobj=open(FILESDICT[db],'rb') ; 
            ATTRDICT[db] =  pickle.load(fobj) ;
            fobj.close() ; 
        elif ATTRDICT[db] :
            pass ; 
        else : 
            raise ValueError('Argument is not a module attribute (mouse/human) and is not a file path.') ; 
            
    except KeyboardInterrupt : 
        print('Load Interrupted.') ; 

    refroots() ; 

def reload(db) : 

    fobj=open(FILESDICT[db],'rb') ; 
    ATTRDICT[db] =  pickle.load(fobj) ;
    fobj.close() ; 

    refroots() ; 

def loadall() : 

    for a in ATTRDICT : 
        load(a) ; 

    refroots() ; 

def rapid(theDict,outfields=('EID','Symbol','Taxon')) : 

    if theDict == None : 
        sys.stderr.write('Reference dictionary not initialized.\n') ; 
        return ; 


    KEYS=theDict.keys()

    def complete(text,state) : 
        for cmd in KEYS : 
            if cmd.startswith(text) : 
                if not state : 
                    return cmd ; 
                else : 
                    state -= 1 ;

    first=True ;
    line='' 

    readline.parse_and_bind("tab: complete") ; 
    readline.set_completer(complete) ; 

    while first or line.strip() : 

        try : 

            first=False ; 

            line=input('==> ') ;
            linel=line.split()

            by=theDict.get(linel[0]) ; 
            if by :
                out=by.get(linel[1]) ; 
            else : 
                continue ; 

            if not out : 
                continue ; 

            outstr=OUTCOLOR+OUTSTYLE

            if type(out) is not list : 
                out=[out] ; 
            if type(out) is None : 
                continue ; 

            for o in out : 
                for f in outfields : 
                    if not o.get(f): 
                        outstr += 'N/A' ; 
                    else : 
                        outstr += '{: <12}'.format(o[f]) ; 
                outstr += '\n' ; 


            sys.stdout.write(outstr+colorama.Style.RESET_ALL) ; 


        except IndexError : 
            break ; 

def update_dup( newdup ) :

    outf=open(FILESDICT['dup'],'wb') ;

    pickle.dump(newdup,outf) ; 

    outf.close() ; 

def load_biogrid(filters='default',m2h=False,h2m=False,force_qualify=True,superdebug=False,qualify='') : 

    global biogrid ; 

    bgpath= REFERENCEPATH + 'biogrid_latest' ; 

    from lib import interactors as I 
    from lib import interactors_extras as ie

    mynf=None ; 
    myif=None ; 
    (myif,mynf)=ie.bg_regex_assembler() ; 

    if filters == 'default' : 
        pass ; 
    elif filters.lower() in 'exogenous' : 
        mynf=ie.exogenous_regex_assembler() ; 
    elif filters.lower() in 'ionly' : 
        mynf=None ; 
    elif filters.lower() in 'physical' : 
        physical_subf=ie.subfilter(r'genetic','systemType') ;
        myif=ie.Filter(excludes={physical_subf,}) ; 
        mynf=None ; 
    elif filters is None or filters.lower() in 'none' : 
        mynf=None ; 
        myif=None ; 
        
    biogrid=I.dataSet(i_filter=myif,n_filter=mynf,superdebug=superdebug) ; 
    bgfile=open(bgpath) ; 
    biogrid.parse(bgfile,fd=I.fd_biogrid,m2h=m2h,h2m=h2m,force_qualify=force_qualify,qualify=qualify) ;
    bgfile.close()

def load_preppi() : 
    global preppi ; 
    from lib import interactors as I 
    from lib import interactors_extras as ie
    preppipath= REFERENCEPATH + 'preppi_150727_lr600.i' ; 
    preppi=I.dataSet() ; 
    preppi_f=open(preppipath) ;
    preppi.parse(preppi_f,fd=I.fdms) ; 
    preppi_f.close() ; 

exores={\
r"HRNR.*",\
r"KRT[0-9]+.*",\
r"Ker.*",\
r"col[0-9]+[AB][0-9]+",\
r'dcd.*',\
r'dsp.*',\
r'flg.*',\
r'DSG.*',\
r'JUP.*',\
r'CSN2',\
r'GFP.*',\
r'HLA.*',\
r'IG[HLEK].*',\
r'ALB',\
r'LYZ',\
r'DSC.*',\
r'CRAP.*',\
r'KRT.*',\
r'CONTAM.*',\
}

endores={\
r"HSP.*",\
r"UB[CQLB].*",\
r"UBA52",\
r"SUMO",\
r"NPM.*",\
r"SPT.N.*",\
r"dnaj.*",\
r'CAD.*',\
r'RPS27A.*',\
r'..RS',\
r'BAG2',\
r'E[EIT]F.*',\
r'RP[SL].*',\
} ;
 
cytores={\
r"TUB[ABG].*",\
r"ACT[BG].*",\
r"MY[OH].*",\
r'FLN[ABC].*',\
r'VIM.*',\
r'ruvbl.*',\
}

highly_variable={\
r'HIST.*',\
r'CSNK.*',\
}

pseudogenes={\
r'LOC.*',\
r'Gm[0-9]\+',\
}

exogenous_regexes=exores

def firstnp(sym,ghash) : 

    for pepacc in hsg['Symbol'][sym]['Peptide'] : 

        if pepacc[0:3] == 'NP_' : 
            return pepacc ;

    else : 
        return None ; 

default_regexes= exores | endores | cytores | highly_variable ; 

def load_refsuite(filters='default',m2h=False,h2m=False,superdebug=False,force=False) : 

    from lib import interactors as I
    from lib import interactors_extras as ie

    global refsuite ;

    if refsuite is not None and len(refsuite.nodes) > 0 and not force : 
        return ; 

    bgpath  =   REFERENCEPATH + 'biogrid_latest' ; 
    empath  =   REFERENCEPATH + '/complexes/emiliome.i' ; 
    bppath  =   REFERENCEPATH + '/bioplex.i' ; 

    mynf=None ; 
    myif=None ; 
    (myif,mynf)=ie.bg_regex_assembler() ; 

    if filters == 'default' : 
        pass ; 
    elif filters.lower() in 'exogenous' : 
        mynf=ie.exogenous_regex_assembler() ; 
    elif filters.lower() in 'ionly' : 
        mynf=None ; 
    elif filters is None or filters.lower() in 'none' : 
        mynf=None ; 
        myif=None ; 

    refsuite=I.dataSet(i_filter=myif,n_filter=mynf,superdebug=superdebug) ; 

    bgfile=open(bgpath) ; 
    refsuite.parse(bgfile,fd=I.fd_biogrid,m2h=m2h,h2m=h2m,force_qualify=True,\
    force_score=0.0,qualify='bg') ; 
    bgfile.close() ; 

    emfile=open(empath) ; 
    refsuite.parse(emfile,fd=I.fd_emili,m2h=m2h,h2m=h2m,force_qualify=True,\
    force_score=0.0,qualify='em') ; 
    emfile.close() ; 

    bpfile=open(bppath) ; 
    refsuite.parse(bpfile,fd=I.fd_emili,m2h=m2h,h2m=h2m,force_qualify=True,directed=True,\
    force_score=0.0,qualify='bp') ; 
    bpfile.close() ; 
    

    #return refsuite.clone() ; 
def load_dompairs() : 

    global dompairs ; 
    dompairs=list() ; 

    dompairsfile= open( REFERENCEPATH + 'bioplex_predicted.tsv' ,'r' ) 

    for line in dompairsfile :  
        linel=line.strip().split('\t') ; 
        dompairs.append(set(linel)) ; 

    dompairsfile.close() ; 

   #def quickconvert(eid) : 

   #    load('m2h') ; 
   #    load('h2m') ; 
   #    load('hmg') ; 

   #    if eid in m2h : 
   #        for tryeid in m2h.get(eid,set()) : 
   #            if rbas

   #    elif eid in h2m : 

   #    else : 
   #        raise IndexError('Supplied eid '+repr(EID)+' could not be found in either conversion dictionary.') ; 

def list_domains(rdict,key,outfields=('Pssm','Name'),root_families=True) : 

    load('cdd') ; 

    global VALID_CDDS ;

    if not VALID_CDDS : 
        VALID_CDDS=set(cdd['ByPssm'].keys()) ; 

    try :
        int(key) ; 
        searchkey='EID' ; 
    except ValueError : 
        searchkey='Symbol' ; 

    if root_families : 
        for c in set(rdict[searchkey][key]['CDD']) & VALID_CDDS : 

            sys.stdout.write('\t'.join([cdd['ByPssm'][c]['Root'][of] for of in outfields])) ; 
            sys.stdout.write('\n') ; 
    else : 
        for c in set(rdict[searchkey][key]['CDD']) & VALID_CDDS : 

            sys.stdout.write('\t'.join([cdd['ByPssm'][c][of] for of in outfields])) ; 
            sys.stdout.write('\n') ; 

def syms_with_domain(rdict,pssm,root_families=True,outfields=(['EID','Symbol'])) : 
    
    load('cdd') ; 

    global VALID_CDDS ;

    if not VALID_CDDS : 
        VALID_CDDS=set(cdd['ByPssm'].keys()) ; 

    if type(pssm) is str : 
        pssm={pssm} ; 

    if root_families :
        for e in rdict['EID'] : 
            if rdict['EID'][e]['CDD'] is not None and \
             any([ cdd['ByPssm'][c]['Root']['Pssm'] in pssm for c in set(rdict['EID'][e]['CDD']) & VALID_CDDS ]) : 
                sys.stdout.write('\t'.join([ rdict['EID'][e][of] for of in outfields ]))
                sys.stdout.write('\n') ; 
    else :
        for e in rdict['EID'] : 
            if rdict['EID'][e]['CDD'] is not None and \
             any([ cdd['ByPssm'][c]['Pssm'] in pssm for c in set(rdict['EID'][e]['CDD']) & VALID_CDDS ]) : 
                sys.stdout.write('\t'.join([ rdict['EID'][e][of] for of in outfields ]))
                sys.stdout.write('\n') ; 
def paperlimit(rdict,limit) : 

    papercounts=dict() ;

    for e in rdict['EID'].keys() : 
        for p in rdict['EID'][e]['Pubmed'] : 
            if p in papercounts : 
                papercounts[p] +=1 ; 
            else :
                papercounts.update({ p : 1 }) ;

    return { p for p in papercounts.keys() if papercounts[p] > limit }
                


def merge(*args) : 

   #def isiterable(obj) : 
   #    try : 
   #        iter(obj)
   #        return True ; 
   #    except TypeError :
   #        return False ; 

   #def recursive_iter(obj,wdict) : 

   #    if type(obj) is not str and isiterable(obj) : 

   #        for o in obj :

   #        return recursive_iter(obj) :

    #1) find all terminal entries
    #2) make new dict with re-indexed entries
    from lib.rbase import VALID_INDEXING_FIELDS 
    terminal_entries=dict() ; 

    for d in args : 
      for k1 in d.keys() :
        for k2 in d[k1].keys() : 
            # e.g. entrez, sym, org

            if type(d[k1][k2]) is list : 
                for l in d[k1][k2] : 
                    terminal_entries.update({l['Entrez']+l['Symbol'] : l}) ; 
            elif type(d[k1][k2]) is dict : 
                terminal_entries.update({ d[k1][k2]['Entrez']+d[k1][k2]['Symbol'] : l }) ; 

    output=dict() ; 
    for v in VALID_INDEXING_FIELDS : 
        output.update({ v : dict() }) ;

    for t in terminal_entries : 
        for vif in VALID_INDEXING_FIELDS : 
            if t.get(vif) : 
                if not output.get(vif) : 
                    output[vif].update({ t.get(vif) : t }) ; 
                elif output.get(t[vif]) and type(output.get(t[vif])) is list : 
                    output[vif][t[vif]].append(t) ; 
                else : 
                    output[vif].update({ t[vif] : [ output[vif][t[vif]] , t ]})  ; 

    return output

def domroots(domlist) : 
    outset=set()
    load('cdd') ; 

    for d in domlist : 
        if d in cdd['ByPssm'] :
            outset.add( cdd['ByPssm'][d]['Root']['Pssm'] ) ;

    return outset

quicksearch=lambda x : rapid(hsg,outfields=('Symbol','Summary')) ; 
def pmid_pairs(id1,id2,ref=hsg) :

    from engarde import BIG_SCREENS
    id1type='s' ; 
    id2type='s' ; 
    try :
        foo=int(id1type)
        id1type='e' ; 
    except ValueError : 
        pass ; 

    try :
        foo=int(id2type)
        id2type='e' ; 
    except ValueError : 
        pass ; 

    if id1type == 's' :   
        pmids1=rbase.hsg['Symbol'][s]['Pubmed'] ; 
    else : 
        pmids1=rbase.hsg['EID'][e]['Pubmed'] ; 

    if id2type == 's' :    
        pmids2=rbase.hsg['Symbol'][s]['Pubmed'] ; 
    else : 
        pmids2=rbase.hsg['EID'][e]['Pubmed'] ; 

    return pmids1  & pmids2 - E.BIG_SCREENS

def pubdump(pmids,**kwargs) : 
    from engarde import pubdump
    pubdump(pmids,**kwargs) ; 

# RESUME : update dict to map terminal entries

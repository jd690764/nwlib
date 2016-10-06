import re
import io
import urllib.request, urllib.error, urllib.parse  
import sys

from Bio import Entrez
from lib import config as cf

#
#   Conventions : 'find' indicates an internet search
#               'fetch' indicates retrieval of information from the internet (e.g given an accession)
#               'get' indicates a record search
#
#               GR: gene record
#               PR: protein record
#               Pub : pubmed record
#
#               Sym : gene symbol
#               EID : entrez gene ID
#               PA : protein accession (NCBI)
#               Sw : protein accession (swissprot/trembl/uniprotkb) ; 
#               Org: NCBI taxonomy ID
#               Sum: NCBI Entrez summary string
#               Seq: sequence (AA or NT)
#               Loc: chromosome and arm (e.g. 12p)
#               Abs : abstract (pubmed record)
#

            # mouse below  the break
BIG_SCREENS={ '25416956','12477932','15489334','14702039','17081983','10470851',\
              '16344560','16189514','21900206','16169070','20360068','19615732',\
              '22939629','11347906','23251661','12853948','15146197','11230166',
              '17353931','23376485','22354994','12690205','23602568',\
              '19056867','9074931','24248522','18251583','21630459',\
              '23580065',\
              '10349636','11076861','16141072','12466851','16141073','12477932',\
              '11042159', '11217851', '15489334','16602821','10922068','21677750',\
              '14681479','23788249','18287559','21905163','19447967',\
              '21988832','15231748','16381901','19199708','11076863','15489336',\
              '12665801','26344197','18029348','10718198','26186194','18854154','22681889',\
              '21516116','22190034','11042154'
               }


OTIDICT={   'Varicella-zoster virus (strain Dumas)' : '10335' ,\
            'Varicella-zoster virus (strain Oka vaccine)' : '10335' ,\
            'Homo sapiens'                          : '9606',\
            'Mus musculus'                          : '10090',\
            'Sus scrofa'                            : '9823',\
            'Rattus norvegicus'                     : '10116',\
            'Bos taurus'                            : '9913',\
            'Gallus gallus'                         : '9031',\
            'Xenopus laevis'                        : '8355',\
            'Schizosaccharomyces pombe'             : '4896',\
            'Saccharomyces cerevisiae'              : '4932',\
            'Escherichia coli'                      : '562',\
            'Caenorhabditis elegans'                : '6239',\
            'Chlamydomonas reinhardtii'             : '3055',\
            'Drosophila melanogaster'               : '7227' }

Entrez.email="mrkelly@stanford.edu"

refToCurrent=re.compile(".*Gene-track_current.*") ; 

def findEID(query,txid="",trySym=True) :


    outl=[] ;

    for symCriterion in ["[sym]",""] :
    # getting ids is tried twice, once with the sym criterion and once without

        if (symCriterion == "[sym]" and not trySym ) :
            continue ; 

        tries=0 ;

        if ( txid ) :
            searchstr=query+symCriterion+" AND "+txid+"[Taxonomy ID]" ;
        else :
            searchstr=query ;

        while tries < 3 :
            try : 
                r=Entrez.esearch(db="gene",term=searchstr) ;
                break ; 
            except IOError : 
                tries += 1 ; 
        if tries == 3 :
            sys.stderr.write("3 strikes!\n") ; 
            raise ;

        idGetter=re.compile(r"<Id>([0-9]*)</Id>") ;

        s=r.readline() ;
        while s :
            m=idGetter.match(s)
            if  m :
                outl.append(m.group(1)) ; 
            s=r.readline() ;

        if ( len(outl) == 1 ) :
            return outl[0] ;
        elif len(outl) == 0  :
            continue ; 
        else : 
            return outl ;

    return None ; 


def findPA(query,txid="") :

    outl=[] ; 

    tries=0 ;

    if ( txid ) :
        osName=list(OTIDICT.keys())[list(OTIDICT.values()).index(txid)] ;
        searchstr=query+" AND "+osName+"[organism]" ;
    else :
        searchstr=query ;

    while tries < 3 :
        try : 
            r=Entrez.esearch(db="protein",term=searchstr,sort="relevance") ;
            break ; 
        except IOError : 
            tries += 1 ; 
    if tries == 3 :
        sys.stderr.write("3 strikes!\n") ; 
        raise ;

    idGetter=re.compile(r"<Id>([0-9]*)</Id>") ;

    s=r.readline() ;
    while s :
        m=idGetter.match(s)
        if  m :
            outl.append(m.group(1)) ; 
        s=r.readline() ;

    if ( len(outl) == 1 ) :
        return outl[0] ;
    elif len(outl) == 0  :
        return None ; 
    else : 
        return outl ;

    return None ; 

#TODO : 
#findTA,
#fetchTA
#TA get (etc.)

def findPub(query,retmax=20) : 
    outs=set() ; 
    tries=0 ;
    while tries < 3 :
        try : 
            r=Entrez.esearch(db="pubmed",term=query,sort="relevance",retmax=retmax) ;
            break ; 
        except IOError : 
            tries += 1 ; 
    if tries == 3 :
        sys.stderr.write("3 strikes!\n") ; 
        raise ;

    idGetter=re.compile(r"<Id>([0-9]*)</Id>") ;

    s=r.readline() ;
    while s :
        m=idGetter.match(s)
        if  m :
            outs.add(m.group(1)) ; 
        s=r.readline() ;

    return outs ;

def fetchPR(qid,currentOnly=True) :

    records=list() ;
    if ( type(qid) not in [list,tuple,set] ) :
        ql=[qid] ;
    else :
        ql=qid ; 
    
    for q in ql : 
    
        # try 3 times to connect to NCBI
        tries=0 ;
        while tries < 3 :
            try :
                fetched=Entrez.efetch(db="protein",id=q,rettype='gp',retmode="text") ; 
                break ; 
            except IOError : 
                # a subclass of this error is raised on connection failures
                tries +=1 ;
        if tries == 3 :
            sys.stderr.write("3 strikes!\n") ; 
            raise ;
    
        r=io.StringIO() ;
        # prepare a stringIO object (rewindable, unlike einfourl)
        s=fetched.readline() ;
        while s : 
            r.write(s) ; 
            s=fetched.readline() ;
    
        if ( currentOnly ) :
    
            recordIsCurrent=True ;
            r.seek(0) ;
    
            for line in r :
                if ( refToCurrent.match(line) ) :
                    recordIsCurrent=False ;
                    break ;
    
            r.seek(0)
    
            if recordIsCurrent :
                records.append(r) ; 
    
        else : 
            records.append(r) ;
    
    
    if len(records)==0 :
        return ;
    elif len(records)==1 :
        return records[0] ;
    else:
        return records ; 

def fetchGR(qid,currentOnly=True) :

    records=list() ;

    if ( type(qid) not in [list,tuple,set] ) :
        ql=[qid] ;
    else :
        ql=qid ; 

    for q in ql : 

        # try 3 times to connect to NCBI
        tries=0 ;
        while tries < 3 :
            try :
                fetched=Entrez.efetch(db="gene",id=q,retmode="xml") ; 
                break ; 
            except IOError : 
                # a subclass of this error is raised on connection failures
                tries +=1 ;
        if tries == 3 :
            sys.stderr.write("3 strikes!\n") ; 
            raise ;

        r=io.StringIO() ;
        # prepare a stringIO object (rewindable, unlike einfourl)
        s=fetched.readline() ;
        while s : 
            r.write(s) ; 
            s=fetched.readline() ;

        # 

        if ( currentOnly ) :

            recordIsCurrent=True ;
            r.seek(0) ;

            for line in r :
                if ( refToCurrent.match(line) ) :
                    recordIsCurrent=False ;
                    break ;

            r.seek(0)

            if recordIsCurrent :
                records.append(r) ; 

        else : 
            records.append(r) ;


    if len(records)==0 :
        return ;
    elif len(records)==1 :
        return records[0] ;
    else:
        return records ; 

def GRgetSym(record) :


    def symfromrecord(r) : 
        gs=re.compile(r".*Gene-ref_locus>([^<]*)<.*") ;

        r.seek(0) ;
        s=r.readline() 
        while s : 
            m=gs.match(s) ; 
            if ( m ) :
                return m.group(1) ;
            s=r.readline() ;

        return None ;

    if type(record) in [list,tuple,set]  : 
        out=list() ;
        for r in record :
            out.append(symfromrecord(r)) ;
        return out ; 
    else :
        out=symfromrecord(record) ;
        return out ;

def GRgetEID(record) : 

    def eidfromrecord(r) :
        geid=re.compile(r".*Gene-track_geneid>([^<]*)<.*") ;

        r.seek(0) ;
        s=r.readline() ;

        while s : 
            m=geid.match(s) ;
            if ( m) :
                return m.group(1) ; 
            s=r.readline() ;

        return None ;

    if type(record) in [list,tuple,set]  : 
        out=list() ;
        for r in record :
            out.append(eidfromrecord(r)) ;
        return out ; 

    else :
        out=eidfromrecord(record) ;
        return out ;

def PRgetEID(protrecord) : 

    def fromrecord(r) :
        eidfetcher=re.compile(r'.*\/db_xref="GeneID:([0-9]*).*"') ; 

        r.seek(0) ;
        s=r.readline() ;

        while s : 
            m=eidfetcher.match(s) ;
            if ( m) :
                return m.group(1) ; 
            s=r.readline() ;

        return None ;

    if type(protrecord) in [list,tuple,set]  : 
        out=list() ;
        for r in protrecord :
            out.append(fromrecord(r)) ;
        return out ; 

    else :
        out=fromrecord(protrecord) ;
        return out ;

def PRgetSym(protrecord) : 

    def fromrecord(r) :
        gsfetcher=re.compile(r'.*\/gene="(.*)"') ; 

        r.seek(0) ;
        s=r.readline() ;

        while s : 
            m=gsfetcher.match(s) ;
            if ( m) :
                return m.group(1) ; 
            s=r.readline() ;

        return None ;

    if type(protrecord) in [list,tuple,set]  : 
        out=list() ;
        for r in protrecord :
            r.seek(0)
            out.append(fromrecord(r)) ;
        return out ; 

    else :
        protrecord.seek(0)
        out=fromrecord(protrecord) ;
        return out ;

def PRgetOrg(protrecord) : 
    def fromrecord(r) :
        orgfetcher=re.compile(r'.*ORGANISM *([A-Z].*)$') ; 

        r.seek(0) ;
        s=r.readline() ;

        while s : 
            m=orgfetcher.match(s) ;
            if ( m) :
                return m.group(1) ; 
            s=r.readline() ;

        return None ;

    if type(protrecord) in [list,tuple,set]  : 
        out=list() ;
        for r in protrecord :
            r.seek(0)
            out.append(fromrecord(r)) ;
        return OTIDICT[out] ; 

    else :
        protrecord.seek(0)
        out=fromrecord(protrecord) ;
        return OTIDICT[out] ;


def GRgetPA(record,productType="peptide",crapFilter=True) :

    # matches the gene commentary type field ( because associated genomic and rna sequences also have accessions)
    # matches accession fields

    def getgp(record,productType,crapFilter) : 
    #get Gene-commentary entries 
        out_acc=set() ;
        gct=re.compile(r".*Gene-commentary_type value=\"([a-z]*)\">.*") ;
        gca=re.compile(r".*Gene-commentary_accession>([A-Za-z0-9_-]*)<.*") ;
        if ( crapFilter) :
            notcrap=re.compile(r".*NP_.*",re.IGNORECASE) ; 

        record.seek(0) ;
        s=record.readline() ;
        while s :

            m_type=gct.match(s) ;
            # identify a  peptide-type gene commentary
            if ( m_type and m_type.group(1) == productType ) :
                #print "\"peptide\" found."
                s=record.readline() ;
                m_acc=gca.match(s) ; 
                m_type=gct.match(s) ;
                # search for the next accession OR gene commentary type

                while s and ( not m_acc  and not m_type)   :
                    s=record.readline() ;
                    #print s ;
                    m_acc=gca.match(s) ;
                    m_type=gct.match(s) ;

                if ( m_acc ) :
                    #print "Found accession "+m_acc.group(1) ;
                    if not crapFilter :  
                        out_acc.add(m_acc.group(1)) ;
                    else : 
                        m=notcrap.match(m_acc.group(1)) ;
                        if m : 
                            out_acc.add(m.group(0)) ;
                elif ( m_type ): 
                    #print "WARNING: passing (type) "+m_type.group(0)   ; 
                    pass ;
                else :
                    #print "WARNING: passing" 
                    pass ;
                # append that accession to the retval

            s=record.readline() ;
            # continue searching until EOF

        return out_acc ; 

    if ( type(record) in [list,tuple,set] ) :
        out=list() ;
        for r in record : 
            out.append(getgp(r,productType,crapFilter)) ; 
        return out ; 
    else :
        return getgp(record,productType,crapFilter) ; 



def closeRecord(record) :

    if ( type(record) in [list,tuple,set] ):
        for r in record:
            r.close() ;
    else : 
        record.close() ;

def GRgetSw(record) :

    upkbHeading=re.compile(r".*<Dbtag_db>UniProtKB/Swiss-Prot.*",re.IGNORECASE) ; 
    goods=re.compile(r".*<Object-id_str>([^>]*)</Object-id_str>.*",re.IGNORECASE) ; 

    record.seek(0) ; 

    scanForGoods=False ;

    for line in record : 

        if not scanForGoods : 
            if upkbHeading.match(line) :
                scanForGoods=True ; 

        else : 
            m=goods.match(line) ; 
            if m : 
                return m.group(1) ; 

    return ;

def GRgetOrg(record) :

    txHeading=re.compile(r".*<Dbtag_db>taxon.*",re.IGNORECASE) ; 
    goods=re.compile(r".*<Object-id_id>([^<]*)</Object-id_id>.*",re.IGNORECASE) ; 

    record.seek(0) ; 

    scanForGoods=False ;

    for line in record : 

        if not scanForGoods : 
            if txHeading.match(line) :
                scanForGoods=True ; 

        else : 
            m=goods.match(line) ; 
            if m : 
                return m.group(1) ; 

    return ;

def fetchSw(swiss,asTuple=False,timeout=5) :
    """
    Swisseq returns a single-line header string by default, or a full FASTA entry as a
    cStringIO object when called with fullSequence=True  

    This function runs INDEPENDENTLY of any Entrez records or methods.

    """

    tries=0 ;

    while tries < 3 :
        try : 
            ur=urllib.request.urlopen(cf.uniprotUrl+swiss+'.fasta',timeout=timeout) ;
            break ; 
        except IOError : 
            tries += 1; 
    if tries == 3 : 
        sys.stderr.write("3 strikes!\n") ; 
        raise ; 


    if not asTuple :
        return ur.readline().strip() ;
    else:
        #r=cStringIO.StringIO() ;
        acc=ur.readline().decode().strip() ; 
        seq=str() ; 
        for line in ur :
            seq += line.decode().strip() ; 
        #return r ; 
        return (acc,seq) ; 

#def recordTimedOut(record) : 


def GRgetSum(records) :

    def fromrecord(r) :
        sumfetcher=re.compile(r'.*Entrezgene_summary>([^<]*)<.*') ; 

        r.seek(0) ;
        s=r.readline() ;

        while s : 
            m=sumfetcher.match(s) ;
            if ( m) :
                return m.group(1) ; 
            s=r.readline() ;

        return None ;

    if type(records) in [list,tuple,set]  : 
        out=list() ;
        for r in records :
            r.seek(0)
            out.append(r) ;
        return out ; 

    else :
        records.seek(0)
        out=fromrecord(records) ;
        return out ;

def PAfetchSeq(recids) :

    def getseq(recid) :
        trials=0 

        while trials < 3 : 
            try :
                record=Entrez.efetch(db='protein',id=recid,rettype='fasta',retmode='text') ;
                break ; 
            except IOError : 
                trials += 1 ;

        if trials==3 :
            raise ; 

        outseq="" ;
        record.readline() ; # gets past the header
        s=record.readline() ;

        while s :
            outseq += s.strip()  ; 
            s=record.readline() ;

        return outseq ; 


    if type(recids) in (list,tuple,set) : 

        outl=[] ;

        for r in recids : 

            outl.append(getseq(r)) ;

        return outl ; 

    else : 
        return getseq(recids) ;

def GRgetLoc(record) : 

    def locfromrecord(r) :

        gl=re.compile(r".*Gene-ref_maploc>([^<]*)<.*") ; 

        r.seek(0) ;
        s=r.readline() ; 
        while s: 
            m=gl.match(s) ;
            if (m) :
                return m.group(1) ; 
            s=r.readline() ; 

        return None


    if type(record) in [list,tuple,set] :
        out=list() ;
        for r in record : 
            out.append(locfromrecord(r)) ; 
        return out ; 
    else:
        out=locfromrecord(record) ; 
        return out ; 

def r2sio(record)  : 

    out=io.StringIO() ;

    for line in record : 
        out.write(line) ;
    record.close() ;

    out.seek(0) ;

    return out ; 

def fetchPub(qid) : 
    records=list() ; 

    if type(qid) not in [list,tuple,set]  : 
        ql=[qid] ; 
    else : 
        ql=qid ; 

    for q in ql : 
        tries=0 ;
        while tries < 3 :
            try :
                fetched=Entrez.efetch(db="pubmed",id=str(q),retmode="xml") ; 
                break ; 
            except IOError : 
                # a subclass of this error is raised on connection failures
                tries +=1 ;
        if tries == 3 :
            sys.stderr.write("3 strikes!\n") ; 
            raise ;
    
        r=io.StringIO() ;
        # prepare a stringIO object (rewindable, unlike einfourl)
        s=fetched.readline() ;
        while s : 
            r.write(s) ; 
            s=fetched.readline() ;
    
        records.append(r) ;
    
    
    if len(records)==0 :
        return ;
    elif len(records)==1 :
        return records[0] ;
    else:
        return records ; 
        

def pubGetAbs(rec) : 

    from lxml import etree
    from textwrap import fill

    rec.seek(0)
    thetree=etree.parse(rec)
    abstract_items=thetree.findall('.//AbstractText')

    if not abstract_items : 
        return 'N/A' ; 

    return str(''.join([a.text for a in abstract_items])) ; 

    
def pubGetTitle(rec) : 

    from lxml import etree
    from textwrap import fill

    rec.seek(0)
    thetree=etree.parse(rec)
    title_items=thetree.findall('.//ArticleTitle')

    if not title_items : 
        return 'N/A' ; 

    return str(''.join([a.text for a in title_items])) ; 

def pubdump(pmids,of=sys.stdout,screens=False) : 
    ws=set(pmids)
    if not screens : 
        ws -= BIG_SCREENS

    for p in ws : 
        rec=fetchPub(p) ; 
        of.write('{:+<60}\n'.format('')) ;
        of.write(p+':  '+pubGetTitle(rec)+'\n')
        of.write('{:~<60}\n'.format('')) ;
        of.write(pubGetAbs(rec)+'\n\n\n')
        of.write('{:+<60}\n'.format('')) ;

def pubGetDate(rec) : 
    from lxml import etree

    rec.seek(0) ; 
    thetree=etree.parse(rec)

    if thetree.findall('.//PubmedBookArticle') : 
        # for books
        year=thetree.findall('.//PubDate')[0].text
        month='12'
        day='31'
    else : 
        datetree=thetree.findall('.//DateCreated')[0] ; 
        year,month,day=[ datetree[i].text for i in range(len(datetree)) ]

    return year+month+day

def publication_history(pmids,of=sys.stdout,screens=False) : 

    #allrecs=list() ;
    from textwrap import fill

    recs=fetchPubBatch(pmids,screens=screens) ;

    rectups=[ (rec,pubGetDate(rec),pubGetPMID(rec)) for rec in recs ] ; 

    srecs=sorted(rectups,key=lambda x : x[1]) ;

    for sr in srecs : 
        fulltitle   =  pubGetTitle(sr[0])

        hrdate=sr[1][0:4]+'-'+sr[1][4:6]+'-'+sr[1][6:] 

        of.write(hrdate+' : '+sr[2]+'\n'+fill(fulltitle)+'\n'  ) ; 

def pubGetPMID(rec) : 
    from lxml import etree
    rec.seek(0) 
    thetree=etree.parse(rec)
    return thetree.findall('.//PMID')[0].text


def fetchPubBatch(pmids,screens=False) : 
    ws=set(pmids)
    if not screens : 
        ws -= BIG_SCREENS

    from multiprocessing import Pool
    thepool=Pool(processes=24)
    recs=thepool.map(fetchPub,ws)
    thepool.close()
    thepool.join()


    for r in recs : 
        r.seek(0) ; 
    return recs
    


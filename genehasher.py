from lxml import etree 
import pickle  
import sys 
from os.path import isfile
import urllib.request, urllib.error, urllib.parse
import io 
import gzip
import re

from lib import markutils as mu
from lib import config as cf

VALID_INDEXING_FIELDS={ 'Symbol','EID','External','Synonym','SwissProt','Peptide','TrEMBL', 'mRNA' }
EMPTY_CD_DICT={ 'Super' : None , 'Sub' : None, 'Pssm' : None , 'Acc' : None , 'Name' : None , 'Root' : None , 'Desc' : None} ;
VIF_PEP={ 'Acc','GI','EID','Sym', 'mRNA' } ;

def eseek(element,childTag) : 

    for c in list(element) : 
        if c.tag == childTag : 
            return c ; 
    else : 
        raise KeyError ; 

def hasAcc(element) :

    if 'Gene-commentary_accession' in [x.tag for x in list(element.getparent())] :
        return True
    else:
        return False ; 

isProt  = lambda x : x.get('value') == 'peptide' ;
isRNA   = lambda x : x.get('value') == 'mRNA'
getAcc  = lambda x : eseek(x.getparent(),'Gene-commentary_accession').text ;
gettext = lambda x : x.text ;

def getProtAccs(ele) : 

    #print map(lambda x : x.get('value'),ele)

    accs=list(map(getAcc,list(filter(hasAcc,list(filter(isProt,list(ele)))))))

    #print accs ;

    alreadyIn=set() ; 
    i=0 ; 
    while i < len(accs) : 
        if accs[i] not in alreadyIn :
            alreadyIn.add(accs[i]) ; 
            i += 1; 
        else : 
            accs.pop(i) ; 

    #print accs  ;
    return accs ;

def getRNAAccs(ele) :

    accs      = list(map(getAcc,list(filter(hasAcc,list(filter(isRNA,list(ele)))))))

    alreadyIn = set() ;
    i=0 ;
    while i < len(accs) :
        if accs[i] not in alreadyIn :
            alreadyIn.add(accs[i]) ;
            i += 1;
        else :
            accs.pop(i) ;
            
    return accs ;
                                                                        

class genehash(dict) : 

    def __init__(self,inTerm=True) :
        super(genehash,self).__init__() ;

        for field in VALID_INDEXING_FIELDS : 
            self.update({ field : dict() }) ;

    def parse(self,infilename,taxon='') : 
        # you should be able to parse multiple times with multiple input files.
        # then pickle the class.

        # the 'taxon' argument is there to constrain entries that are unique to a certain
        # taxonomy id rather than species-- this became an issue distinguishing between modern
        # homo sapiens and the denisovans and neanderthals

        infile=open(infilename,'rb') ;

        def stuff(field,newkey,ng) : 

            if not self[field].get(newkey) : 
                self[field].update({ newkey : ng }) ; 
                # add new item to field that does not yet exist
            elif type(self[field].get(newkey)) is list : 
                self[field][newkey].append(ng);
                # append new gene to list already contained in self[field][newkey]
            else : 
                self[field].update({ newkey : [ self[field][newkey], ng ] }) ;
                # init list with previous contents and with new gene

        #scouting
        numgenes=0 ; 
        while True : 
            try : 
                _=pickle.load(infile) ;
                numgenes += 1 ; 
            except EOFError : 
                break ; 

        infile.seek(0) ;
        ongene=0 ; 
        while True : 
            try : 
                ng=pickle.load(infile) ;

                if taxon and ng['Taxon'] != taxon :
                    sys.stderr.write('NOTE: Taxonomy filtering rejects {}:{} from taxon {} (must be {}).\n'\
                     .format(ng['EID'],ng['Symbol'],ng['Taxon'],taxon)) ;
                    sys.stderr.flush() ;
                    continue ; 


                for field in VALID_INDEXING_FIELDS : 
                    newkey=ng[field] ; 
                    if type(newkey) is list : 
                        for k in newkey : 
                            stuff(field,k,ng) ; 
                    else:
                        stuff(field,newkey,ng) ; 

                ongene += 1; 
                mu.waitbar(80*ongene/numgenes,80,showPct=True) ; 
                
            except EOFError : 
                break ; 

        mu.waitbar(80,80,showPct=True) ; 
        infile.close() ; 
################################################################################

class pephash(dict) : 
    def __init__(self,inTerm=True) :
        super(pephash,self).__init__() ;

        for field in VIF_PEP : 
            self.update({ field : dict() }) ;

    def parse(self,infilename) : 
    # you'll probably have to do this several times.

        def stuff(field,newkey,np) : 

            if not self[field].get(newkey) : 
                self[field].update({ newkey : np }) ; 
                # add new item to field that does not yet exist
            elif type(self[field].get(newkey)) is list : 
                self[field][newkey].append(np);
                # append new gene to list already contained in self[field][newkey]
            else : 
                self[field].update({ newkey : [ self[field][newkey], np ] }) ;
                # init list with previous contents and with new gene

        infile=open(infilename,'r') ;
        npep=0 ;
        for line in infile : 
            if line[0:5] == 'LOCUS' :
                npep += 1; 

        infile.seek(0) ;

        onpep=0 ; 
        np=dict() ;
        seqing=False ; 

        acc='' ;
        peplen=0 ; 
        definition='' ; 
        seq='' ; 
        eid='' ; 
        gi='' ; 
        sym='' ; 
        transcript = ''
        
        for line in infile : 

            if line[0:2]=='//' :

                np={    'Acc'       : acc ,\
                        'GI'        : gi ,\
                        'Length'    : int(peplen),\
                        'Seq'       : seq,\
                        'EID'       : eid,\
                        'Sym'       : sym,\
                        'Def'       : definition.strip(),\
                        'CDs'       : cds,
                        'mRNA'      : transcript } ;

                seqing=False ;

                for field in VIF_PEP : 
                    nk=np[field] ;
                    if type(nk) is list : 
                        for k in nk : 
                            stuff(field,k,np) ;
                    else : 
                        stuff(field,nk,np) ;

                onpep += 1 ;

            elif line[0:5] == 'LOCUS' :

                acc='' ;
                peplen=0 ; 
                definition='' ; 
                seq='' ; 
                eid='' ; 
                gi='' ; 
                sym='' ;
                transcript = ''
                cds=set() ;

                linel=' '.join(line.split()).split() ;
                #acc=linel[1] ;
                peplen=linel[2]

            elif line.startswith('ACCESSION') :
                linel=line.split() ;
                acc=linel[1]

            elif line[0:10] == 'DEFINITION' : 

                for c in line[10:len(line)] : 
                    if c != '[' :
                        definition += c ;
                    else:
                        break ; 

            elif line[0:7] == 'VERSION' :

                #linel=' '.join(line.split()).split() ;
                gi=line.split(':')[1].strip() ;

            elif '/gene=' in line : 
                sym=line.split('"')[1] ;

            elif '/db_xref="GeneID' in line : 
                inter=line.split('"')[1] ;
                eid=inter.split(':')[1] ;

            elif '/db_xref="CDD' in line : 
                inter=line.split('"')[1] ;
                cds.add(inter.split(':')[1]) ;

            elif '/coded_by' in line :
                inter = line.split( '"' )[1]
                transcript = inter.split( '.' )[0]

            elif line[0:6] == 'ORIGIN' : 
                seqing=True ;

            elif seqing : 
                for c in line : 
                    if c.isalpha() :
                        seq += c;

            if npep >= onpep : 
                mu.waitbar(80 * onpep / npep,80,showPct=True) ; 

def parse( infilename ) : 
    # parse genbank file to get protein attributes
    
        infile     = open(infilename,'r')
        onpep      = 0 
        np         = dict()
        seqing     = False
        
        records    = list()
        acc        = ''
        peplen     = 0 
        definition = '' 
        seq        = '' 
        eid        = '' 
        gi         = '' 
        sym        = '' 
        cds        = set()
        mrna       = ''
        
        for line in infile : 
            line = line.rstrip( )
            line = line.rstrip( '\,\;\.' )

            if line[0:2]=='//' :
                if isinstance( eid, list ):
                    eid = ','.join( eid )
                if isinstance( sym, list ):
                    sym = ','.join( sym )
                if isinstance( cds, list ) or isinstance( cds, set ):
                    cds = ','.join( cds )
                if isinstance( acc, list ):
                    acc = ','.join( acc )
                if isinstance( definition.strip( '\s\n\t\,\;\.' ), list ):
                    definition = ','.join( definition.strip( '\s\n\t\,\;\.' ) )
                if isinstance( seq, list ):
                    seq = ','.join( seq )
                if isinstance( peplen, list ):
                    peplen = ','.join( peplen )
                if isinstance( mrna, list ):
                    mrna = ','.join( mrna )

                np={    'Acc'       : acc ,
                        'GI'        : gi ,
                        'Length'    : int(peplen),
                        'Seq'       : seq,
                        'EID'       : eid,
                        'Sym'       : sym,
                        'Def'       : definition,
                        'CDs'       : cds,
                        'mRNA'      : mrna}

                records.append( np )                
                seqing = False
                onpep += 1

            elif line[0:5] == 'LOCUS' :

                acc        = ''
                peplen     = 0 
                definition = '' 
                seq        = '' 
                eid        = '' 
                gi         = '' 
                sym        = '' 
                cds        = set()
                mrna       = ''
                
                linel      = ' '.join(line.split()).split()
                peplen     = linel[2]

            elif line.startswith('ACCESSION') :
                linel      = line.split()
                acc        = linel[1]

            elif line[0:10] == 'DEFINITION' : 

                for c in line[12:len(line)] : 
                    if c != '[' :
                        definition += c
                    else:
                        break 

            elif line[0:7] == 'VERSION' :

                gi     = line.split(':')[1].strip()

            elif '/gene=' in line : 
                sym    = line.split('"')[1]

            elif '/coded_by=' in line : 
                inter  = line.split('"')[1]
                mrna   = inter.split( '.' )[0]
                mrna   = re.sub( 'join\(', '', mrna )
                if re.search( '^complement\(', mrna ):
                    mrna = mrna + ')'
                
            elif '/db_xref="GeneID' in line : 
                inter  = line.split('"')[1]
                eid    = inter.split(':')[1]

            elif '/db_xref="CCDS' in line : 
                inter  = line.split('"')[1]
                cds.add(inter.split(':')[1])

            elif line[0:6] == 'ORIGIN' : 
                seqing = True

            elif seqing :
                
                for c in line : 
                    if c.isalpha() :
                        seq += c;

        return( records )

################################################################################
def xml2bin( infilename, outfilename, type='pickle', insist_on_living=True)  : 
    genes=list() ; 
    # list of dicts

    xf=open(infilename,'rb') ;
    #xf=open('short1.xml','r') ; 
    if not isfile(outfilename) : 
        outfile=open(outfilename,'wb') ; 
        outfile.close() ; 

    outfile = ''
    if type == 'pickle':
        outfile = open(outfilename,'ab') ;
    elif type == 'tsv':
        outfile = open( outfilename, 'at' )

        #tree=etree.parse(xf) ; 
    #root=tree.getroot() ;
    count = 0
    for event ,element in etree.iterparse(xf,events=('end',),tag='Entrezgene') :
    #for element in root.iter("Entrezgene") : 
        #element=prelement ; 
        #element=etree.fromstring(etree.tostring(prelement))

        if event != 'end' :
           continue ; 

        try :
            genetype = eseek(element,'Entrezgene_type').get('value') #  == 'protein-coding' and
            if ( not insist_on_living or element[0][0][1].get('value') == 'live' )  : 
                newGene=dict() ; 
                extDBrefs=list() ; 
                newGene.update({'External' : extDBrefs })
                newGene.update({ 'genetype' : genetype });
            else :

                element.clear() ;
                for ancestor in element.xpath('ancestor-or-self::*'):
                    while ancestor.getprevious() is not None:
                        del ancestor.getparent()[0]

                continue ; 
        except KeyError : 
            sys.stderr.write('Unclear type.n') ; 

        newGene.update({ 'EID' : element[0][0][0].text });

        # 3 is the Entrezgene_gene child
        newGene.update({ 'Symbol' : element[3][0][0].text }); 

        try : 
            newGene.update({ 'Location' : element[3][0][2].text }); 
        except IndexError : 
            newGene.update({ 'Location' : None }); 


        newGene.update({ 'Taxon' : gettext(element[2][0].find(".//Object-id_id")) }) ; 
        #<Entrezgene_source> 2
        #    <BioSource> 0 
        #      <BioSource_genome value="genomic">1</BioSource_genome>
        #      <BioSource_origin value="natural">1</BioSource_origin>
        #      <BioSource_org> 2 
        #        <Org-ref> 0 
        #          <Org-ref_taxname>Homo sapiens</Org-ref_taxname>
        #          <Org-ref_common>human</Org-ref_common>
        #          <Org-ref_db> 2 
        #            <Dbtag>  0 
        #              <Dbtag_db>taxon</Dbtag_db>
        #              <Dbtag_tag> 1 
        #                <Object-id> 0 
        #                  <Object-id_id>9606</Object-id_id> 0
        #                </Object-id>
        

        try :
        # gene-ref_db
            for c in list(eseek(element[3][0],'Gene-ref_db')) : 
            # 3-0-3 : Gene-ref_db  (i tried, anyway)
                #newGene['External'].update({ c[0].text : c[1][0][0].text }) ;
                newGene['External'].append( c[0].text+":"+c[1][0][0].text) ;
        except (IndexError,KeyError) :
            newGene.update({ 'External' : None }) ;

        try : 
            syns=list() ; 
            for c in list(element[3][0][4]) : 
                syns.append(c.text) ;

            newGene.update({'Synonym' : syns })
        except IndexError : 
            newGene.update({ 'Synonym' : None }) ; 


        # save this one for later-- we may need to do a sub-iteration
        try : 
            for dbitem in element.iterfind('.//Dbtag_db') : 
                if dbitem.text=='UniProtKB/Swiss-Prot' :
                    if 'SwissProt' not in newGene: # get the first one
                        newGene.update({ 'SwissProt' : dbitem.getparent()[1][0][0].text }) ; 
                    break ; 
            else : 
                newGene.update({ 'SwissProt' : None })
        except IndexError : 
            newGene.update({ 'SwissProt' : None })

        try : 
            for dbitem in element.iterfind('.//Dbtag_db') : 
                if dbitem.text=='UniProtKB/TrEMBL' :
                    newGene.update({ 'TrEMBL' : dbitem.getparent()[1][0][0].text }) ; 
                    break ; 
            else : 
                newGene.update({ 'TrEMBL' : None })
        except IndexError : 
            newGene.update({ 'TrEMBL' : None })

        # changed to list
        try : 
            newGene.update({ 'Pubmed' : [ e.text for e in element.findall('.//PubMedId' ) ] })
        except KeyError : 
            newGene.update({ 'Pubmed' : None }) ; 

        try : 
           #newGene.update({ 'Peptide' : \
           # getProtAccs(eseek(element,'Entrezgene_locus').findall(".//Gene-commentary_type"))  +
           # getProtAccs(eseek(element,'Entrezgene_comments').findall(".//Gene-commentary_type")) }) ;  

            newGene.update({ 'Peptide' : getProtAccs(element.findall(".//Gene-commentary_type")) }) ; 

        except KeyError : 
            newGene.update({ 'Peptide' : None }) ; 
        # 8-0-5 : Entrexgene_comments, Gene-commentary, Gene-commentary_products

        try :
            newGene.update({ 'mRNA'   : getRNAAccs(element.findall(".//Gene-commentary_type")) }) ;
        except KeyEror :
            newGene.update({ 'mRNA' : None }) ;

        try : 
            newGene.update({ 'Summary' : eseek(element,'Entrezgene_summary').text })
        except KeyError : 
            newGene.update({ 'Summary' : None }) ;

        try :
            newGene.update({ 'Desc' : eseek(element[3][0],'Gene-ref_desc').text })
        except KeyError : 
            newGene.update({ 'Desc' : None }) ;

        ### NEW AND EXTREMELY BOSS : CDD DOMAINS !!!!1111ONE!!
        try : 
            CDDentries=set() ; 
            for dbitem in element.iterfind('.//Dbtag_db') : 
                if dbitem.text=='CDD' :
                    # get domain pssm-id
                    #CDDentries.add(dbitem.getparent()[1][0][0].text)
                    # get domain accession
                    cddacc = dbitem.getparent().getparent().getparent()[1].text
                    cddacc = re.sub( r'^([^\:]+)\:.*', r'\1', cddacc )
                    CDDentries.add( cddacc )

            if not CDDentries : 
                newGene.update({ 'CDD' : None }) ; 
            else : 
                newGene.update({ 'CDD' : list(CDDentries) }) ; 
        except IndexError : 
            # not super sure why this would ever happen here
            newGene.update({ 'CDD' : None })
            raise ValueError ; 


        # 4 December 17 : chromosomal start & end positions

        commentaries=[ c for c in element.findall(".//Gene-commentary")
        if c[0].get('value') == 'assembly' and c[2].text == 'Reference' ]

        intertrees=[ c.find(".//Seq-interval") for c in commentaries ]

        newGene.update({ 'Coordinates' : 
            { eseek(c,'Gene-commentary_heading').text : \
                (int(inter[0].text),int(inter[1].text))  for c,inter in zip(commentaries,intertrees) }})

        element.clear() ;
        for ancestor in element.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]

        if type == 'pickle':
            pickle.dump(newGene,outfile,protocol=2) ;
        elif type == 'tsv':
            line = ''
            counter = 0
            for k, v in sorted( newGene.items() ):
                if isinstance(v, list):
                    v = ";".join( v )
                elif v == None:
                    v = ''
                if counter == 0:
                    line = v
                else :
                    line = line + "\t" + v
                counter = counter + 1
            line = re.sub( r'\n', '', line ) 
            outfile.write( line + "\n" )
            
        outfile.flush() ; 

    outfile.close() ; 

def r2sio(record)  : 

    out=io.StringIO() ;

    for line in record : 
        out.write(str(line)) ;
    record.close() ;

    out.seek(0) ;

    return out ; 

def mkCDDtree(debug=False) :

    if debug : 
        sys.stdout.write('Fetching family_superfamily_links... ') ; 
        sys.stdout.flush() ; 
    cdFamMap = r2sio(urllib.request.urlopen(cf.superfamilyFtp)) ; 
    if debug : 
        sys.stdout.write('Done\n.') ; 
        sys.stdout.flush() ; 

    if debug : 
        sys.stdout.write('Fetching CDDID table ... ') ; 
    sio = r2sio(urllib.request.urlopen(cf.cddIdTableFtp)) ; 
    if debug : 
        sys.stdout.write('Decompressing ... ') ; 
        sys.stdout.flush() ; 
    cdTbl = gzip.GzipFile(fileobj=sio) ;
    if debug : 
        sys.stdout.write('Done\n.') ; 
        sys.stdout.flush() ; 

    cdTbl.seek(0) ;
    outdict = createCDDtree( cdFamMap, cdTbl, debug )

    cdTbl.close() ; 
    sio.close() ; 

    return outdict
    
def createCDDtree( cdFamMap, cdTbl, debug = False ):    
    
    def stuff(theDict,theKey,theValue) : 

        if not theDict.get(theKey) : 
            theDict.update({ theKey : theValue }) ; 
            # add new item to field that does not yet exist
        elif type(theDict.get(theKey)) is list : 
            theDict[theKey].append(theValue);
            # append new gene to list already contained in self[field][newkey]
        else : 
            theDict.update({ theKey : [ theDict[theKey], theValue ] }) ;
            # init list with previous contents and with new gene

    outdict={ 'ByAcc' : {} , 'ByPssm' : {} } ;
    
    for line in cdTbl : 

        linel = line.strip().split('\t') ;

        entry = dict(EMPTY_CD_DICT)  ; 
        entry['Pssm']=linel[0] ;
        entry['Acc']=linel[1] ; 
        entry['Name']=linel[2] ; 
        entry['Desc']=linel[3] ; 

        outdict['ByPssm'].update({ linel[0] : entry }) ;
        outdict['ByAcc'].update({ linel[1] : entry }) ;

    for line in cdFamMap : 

        # first, identify all single-step parent/child relations

        linel = line.strip().split('\t') ; 
        if debug:
            print( 'family_relations: ' + '+'.join( linel ))
        
        if linel[1] != linel[3] :
            if linel[1] in outdict['ByPssm'] and linel[3] in outdict['ByPssm']:
                outdict['ByPssm'][linel[1]].update({ 'Super' : outdict['ByPssm'][linel[3]] }) ; 
                stuff(outdict['ByPssm'][linel[3]],'Sub',outdict['ByPssm'][linel[1]]) ; 
            if debug:
                print( linel[3] + ' is parent to ' + linel[1] )
        else :
            if linel[1] in outdict['ByPssm']:
                outdict['ByPssm'][linel[1]].update({ 'Super' : None }) ; 
                if debug:
                    print( 'no parent for ' + linel[1] )

    for domain in list(outdict['ByPssm'].values()) : 

        if domain['Super'] == None :
            domain['Root'] = domain ;
            if debug:
                print( domain['Acc'] + ': root = super' )
        else : 
            nd=domain ; 
            while nd['Super'] != None : 
                nd=nd['Super'] ;
            domain['Root'] = nd ;
            if debug:
                print( domain['Acc'] + ' root is ' + nd['Acc'] )
                
    return outdict ;


#class PubHash(object) : 
#
#    Entrez  =   __import__('Bio.Entrez') ; 
#    mp      =   __import__('multiprocessing') ;
#
#    def __init__() :
#
#        self.Entrez.email('foo@stanford.edu') ;
#        self.pmids=set() ;
#        self.recs=list() ; 
#
#    def parse(ghash) : 
#
#        for eid in ghash.eids.keys() : 
#            self.pmids = ghash['EID'][eid]['Pubmed'] ; 
#            # RESUME

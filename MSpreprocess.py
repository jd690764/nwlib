import sys
import scipy as sp
from scipy import stats 
import numpy as np 
import re
import os
import time
from subprocess import  call
import pickle
import openpyxl as ox

import lib.markutils as mu
import lib.peputils as pep 
import lib.engarde as E 
import lib.rbase as rb

from network.models import Entrez, Ncbiprot

tremblre      = re.compile(r'.*tr\|([^\|]{6})\|.*') 
swisspre      = re.compile(r'.*sp\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})\|.*', re.IGNORECASE ) 
entrezre      = re.compile(r'.*ref\|([\d]+)\|.*') 
symbolre      = re.compile(r'.*GN=([A-Za-z0-9]*).*') 
proteinre     = re.compile(r'.*ref\|([ANYXZ]P_\d+)\.?\d{0,2}\|.*') 
PSEUDO_LENGTH = 375.0
orgs          = { 'hs': 9606, 'mm': 10090 }

def eseek(element,childTag) : 

    for c in list(element) : 
        if c.tag == childTag : 
            return c 
    else : 
        raise KeyError 
#
# IMPORTANT : What type of file should you run this on ? 
#
#
#   Answer 1 : Old bill lane files. These are in nothing like a uniform format, so you need to reduce
#       them down to the following
#   Column 1 : identifier in 'pipe syntax' -- looks for 'sp|(something)|' or 
#   gi|(something) or ref|(something), and then goes to search NCBI and/or EMBL for (something)
#
#   column 2 : total counts.
#   AND THEN YOU NEED TO PUT ID NUMBERS BEFORE EACH ROW
#
#   Answer 2 : New SUMS files.
#
#       Column 1: index (1-
#
#       Column 2: description ( mostly refseq, but we'll parse this as we would the identifiers from the
#           Lane datasets
#
#       Columns 3- : peptide counts by fraction
#       
#       Last column-- total.

rb.load('dup')

hsgREF = dict()
mmgREF = dict()
hspREF = dict()
mmpREF = dict()

gREF   = dict()
pREF   = dict()

def make_gene_reference( org ):

    genes = Entrez.objects.filter( taxid = org ).values( 'eid', 'symbol', 'peptide', 'synonym', 'taxid', 'swissprot', 'trembl')
    ref   = { x['eid']: x for x in genes }
    REF   = { 'eid' : ref }

    for k in [ 'symbol', 'peptide', 'synonym', 'swissprot', 'trembl' ]:
        ref  = dict()
        for r in genes:
            if r[k] and r[k] not in ['', '-']:
                ids = r[k].split(';')
                for id in ids:
                    if id in ref: # the key already exists
                        if not isinstance( ref[id], list): # if the value is not yet a list
                            ref[id] = [ ref[id] ]
                        ref[id].append( r )
                    else:
                        ref[id] = r
            else:
                pass

        REF[ k ] = ref
            
    if org == orgs['hs']:
        global hsgREF
        hsgREF      = REF

    elif org == orgs['mm']:
        global mmgREF
        mmgREF      = REF

    return REF

        
def make_protein_reference( org ):

    prot   = Ncbiprot.objects.all().values( 'gi', 'acc', 'eid', 'protname', 'len', 'symbol', 'taxid' )
    if org == orgs[ 'hs' ]:
        global hspREF
        hspREF = { x['acc']: x for x in prot if x['taxid'] == org }
        return hspREF
    else:
        global mmpREF
        mmpREF = { x['acc']: x for x in prot if x['taxid'] == org }
        return mmpREF
        
        
for org in orgs.values():
    gREF[org] = make_gene_reference( org )
    pREF[org] = make_protein_reference( org )
    

mm2hs = '' #rbase.m2h 
hs2mm = '' #rbase.h2m

#for k in hspREF:
#    print( k )
#    for kk in hspREF[k]:
#        print( kk )
#        print( hspREF[k][kk])
#        break

exores  = [ r'HRNR.*',r'KRT[0-9]+.*',r'KER.*',r'col[0-9]+[AB][0-9]+',r'DCD', r'DSP',r'FLG',r'IG[GLK]',r'CONTAM',r'GFP' ] 
exo     = [re.compile(x,re.IGNORECASE) for x in exores]
idrevs  = re.compile(r".*>rev.*|.*>r-.*",re.IGNORECASE)
getgn   = re.compile(r".* GN=([A-Za-z0-9_.-]*).*")
getos   = re.compile(r".* OS=([A-Za-z0-9 ()-]*) ..=.*") 
contam  = re.compile(r".*contam.*|.*>uc.*",re.IGNORECASE)
wastrd  = re.compile(r".*truncate.*",re.IGNORECASE)
spacc   = re.compile(r".*sp\|(.*)\|.*",re.IGNORECASE)
entrezh = re.compile(r'.*gi\|.*',re.IGNORECASE) 

class MSdatum(object) :

    def __init__( self, idn, desc = "", fxncounts = list(), isReverse = False, score = 0.0 ) :
        self.idn          = idn 
        self.desc         = desc 
        self.fxncounts    = fxncounts 
        self.score        = score
        self.official     = ""
        self.entrez       = ""
        self.organism     = ""
        self.peplen       = 0 
        self.isReverse    = isReverse
        self.totalcounts  = sum(fxncounts)
        
        if fxncounts : 
            self.maxcount = max(fxncounts)
        else : 
            self.maxcount = self.totalcounts 


    def setScore(self,score) :
        self.score        = score

    def setEntrez(self,entrez) : 
        self.entrez       = entrez

    def setOrganism(self,organism) :
        self.organism     = organism

    def setOfficial(self,official) :
        self.official     = official 

    def setLen(self,length) : 
        self.peplen       = length 

    def reassess(self) :
        self.totalcounts  = sum(self.fxncounts) 


class MSdata(object) :

    def __init__(self,name,debug=False):
    # cd and ed are correction and entrez dictionaries

        self.name        = name 
        self.bait        = None 
        self.fpd         = None # "false positive distribution"
        self.fwdata      = list()
        self.rvdata      = list()
        self.debug       = debug
        self.infnames    = [] 
        self.__widn__    = 0 
        self.pseudoscore = 0.0   
        self.background  = dict()
        
        # key : symbol
        # value : counts in blanking data

    def parseLane(self,infobj,unTruncate=False,sep='\t') : 
        # argument is now an input file object
        # dec15/jan15 revision : parse from excel

        f          = infobj
        self.infnames.append(f.name) 

        s          = f.readline().strip().split(sep) # moves past headers
        maxCtDatum = None 
        maxCt      = 0 
        for line in f :
            linel    = line.strip().split(sep)
            thisdesc = linel[1]

            if contam.match(thisdesc) :
                continue

            ###### handle when linel[2] =='NaN' RESUME
            try :
                frac_counts = [ int(linel[2]) ]
            except ValueError : 
                frac_counts = [0] 

            if ( idrevs.match(thisdesc) ) :
                isReverse   = True 
            else :
                isReverse   = False 

            if ( unTruncate and not isReverse and wastrd.match(thisdesc)) :
                sys.stderr.write("WARNING:   Description line\n{}\n seems to have been truncated.\n".\
                 format(thisdesc))

                try :
                    os.remove('.tmp.txt')
                except OSError :
                    pass

                if entrezh.match(thisdesc) : 
                    pass 
                else : 
                    thisdesc = seqRetter(thisdesc)
                    sys.stderr.write("            Calling EMBOSS seqret\n")

            if ( self.debug) :
                sys.stderr.write("DEBUG:   Interaction with ID {} :\n".format(linel[0])) 
                sys.stderr.write("DEBUG:       {}\n".format(thisdesc))
                sys.stderr.write("DEBUG:       {!r}\n".format(frac_counts))
                sys.stderr.write("DEBUG:       Id'd as \"Reverse\": {}.\n".format(isReverse))

            newDatum        = MSdatum( idn = int(linel[0]), desc = thisdesc, fxncounts = frac_counts, isReverse = isReverse )

            if maxCt < newDatum.totalcounts :
                maxCt       = newDatum.totalcounts 
                maxCtDatum  = newDatum 

            if isReverse : 
                self.rvdata.append(newDatum)
            else :
                self.fwdata.append(newDatum)

        self.bait  = maxCtDatum 

        f.close()

    def set_background( self, fname, bestpepdb = 'RPHs', reference = hsgREF ) :

        if bestpepdb == 'RPMm':
            reference = mmgREF
        
        self.background = dict() 

        wb = ox.load_workbook(fname) 
        ws = wb['Proteins'] 
        # to protein worksheet
        # identifiers in col 2, # of spectra in column 6 (0-indexed)

        for row in ws.rows[1:] : 

            entrez, sym, org = desc_interpreter( row[1].value, tryhard = True, bestpepdb = bestpepdb, reference = reference ) 

            if not row[6].value : 
                continue 
            if not self.background.get(sym) : 
                self.background.update({ sym : row[6].value }) 
            else : 
                self.background[sym] += row[6].value

    def parseSUMS( self, infobj, unTruncate = True, sep = '\t' ) :

        f            = infobj
        s            = f.readline().strip().split(sep) # moves past headers
        maxCtDatum   = None 
        maxCt        = 0 
        self.infnames.append(f.name) 
        
        for line in f :
            
            linel    = line.strip().split(sep)
            thisdesc = linel[1]

            if contam.match(thisdesc) :
                continue

            frac_counts = list()

            for c in linel[2:len(linel)-2] :
            # this is MINUS 2 because the second-to-last column is MAX and should not be summed
            # with the others 

                try : 
                    frac_counts.append(int(c))
                except ValueError : 
                    frac_counts.append(0)

            if ( idrevs.match(linel[1]) ) :
                isReverse = True 
            else :
                isReverse = False 

            if ( unTruncate and not isReverse and wastrd.match(thisdesc)) :
                sys.stderr.write("WARNING:   Description line\n{}\n seems to have been truncated.\n".\
                 format(thisdesc))

            if ( self.debug) :
                sys.stderr.write("DEBUG:   Interaction with ID {} :\n".format(linel[0])) 
                sys.stderr.write("DEBUG:       {}\n".format(thisdesc))
                sys.stderr.write("DEBUG:       {!r}\n".format(frac_counts))
                sys.stderr.write("DEBUG:       Id'd as \"Reverse\": {}.\n".format(isReverse))

            newDatum = MSdatum( idn = int(linel[0]), desc = thisdesc, fxncounts = frac_counts, isReverse = isReverse )

            if maxCt < newDatum.totalcounts :
                maxCt      = newDatum.totalcounts 
                maxCtDatum = newDatum 

            if isReverse : 
                self.rvdata.append(newDatum)
            else :
                self.fwdata.append(newDatum)

        self.bait = maxCtDatum 

        f.close()

    def parse_pub(self,sep='\t') : 

        f=infobj
        self.infnames.append(f.name) 

        s=f.readline().strip().split(sep) # moves past headers
        nFields=len(s)

        maxCtDatum=None 
        maxCt=0 

        for line in f :
            linel=line.strip().split(sep)

            thisdesc='> '+linel[3]+' GN='+linel[7]

            if contam.match(thisdesc) :
                continue

            frac_counts=[ int(linel[9]) ];

            if ( idrevs.match(linel[1]) ) :
                isReverse=True 
                if ( self.debug) :
                    sys.stderr.write("DEBUG:   Interaction with ID {} :\n".format(linel[0])) 
                    sys.stderr.write("DEBUG:       {}\n".format(thisdesc))
                    sys.stderr.write("DEBUG:       {!r}\n".format(frac_counts))
                    sys.stderr.write("DEBUG:       Id'd as \"Reverse\": {}.\n".format(isReverse))
            else :
                isReverse=False 


            newDatum=MSdatum(idn=int(linel[0]),desc=thisdesc,fxncounts=frac_counts,isReverse=isReverse)

            if maxCt < newDatum.totalcounts :
                maxCt=newDatum.totalcounts 
                maxCtDatum=newDatum 

            if isReverse : 
                self.rvdata.append(newDatum)
            else :
                self.fwdata.append(newDatum)

        if not self.bait : 
            self.bait=maxCtDatum 

        f.close()


    def parseXML(self,infobj) : 

        from lxml import etree

        xf=infobj
        self.infnames.append(xf.name) 

        if not self.fwdata : 
            data_to_add=dict() 
        else : 
            data_to_add=dict(list(zip( [x.desc for x in self.fwdata] ,self.fwdata))) 
            self.fwdata=list() 

        for event,element in etree.iterparse(xf,events=('end',),tag='spectrum') :

            if event != 'end' : continue 

            thisdesc='' 
            thisref='' 

            try :
                thisdesc=eseek(element,'description').text  
            except KeyError : 
                pass 

            try : 
                thisref=eseek(element,'ref').text                 
            except KeyError : 
                pass 

            if thisdesc or thisref : 

                descline='>'+thisref+' '+thisdesc 

                if descline not in data_to_add : 
                    if idrevs.match(descline) : 
                        isReverse=True 
                    else : 
                        isReverse=False 

                    self.__widn__ += 1 
                    data_to_add.update({ descline : \
                     MSdatum(idn=self.__widn__,desc=descline,fxncounts=[ 1 ],isReverse=isReverse) }) 

                else : 
                    data_to_add[descline].fxncounts[0] += 1 


        maxct=0 
        maxctdatum=None 


        for key,val in list(data_to_add.items()) : 

            if val.isReverse : 

                self.rvdata.append(val) 

            else : 

                val.reassess()

                if val.totalcounts > maxct : 
                    maxct = val.totalcounts 
                    maxctdatum=val 

                self.fwdata.append(val) 

        if not self.bait : 

            self.bait=maxctdatum 


        xf.close() 
                    

    def inferFP(self,scoreBy="max") :
        #thiscomment
        if scoreBy in ( "max","Max","MAX" ) :
            p=lambda x : x.maxcount  -1
        elif scoreBy in ( "total","Total","TOTAL" ):
            p=lambda x : x.totalcounts -1
        else :
            sys.stderr.write("SCOREBYFP: Unsure what to score ")
            return

        inputs=list()

        for d in self.rvdata :
            inputs.append(p(d))

        self.fpd=stats.poisson(np.mean(inputs))


    def score(self,scoreBy="total",norm2len=True) :

        if scoreBy in ( "max","Max","MAX" ) :
            p=lambda x : x.maxcount
        elif scoreBy in ( "total","Total","TOTAL" ):
            p=lambda x : x.totalcounts
        else :
            sys.stderr.write("SCOREBYFP: Unsure what to score ")
            return

        def notNone(serie) : 
            for x in serie : 
                if x != None : 
                    return True
            else : 
                return False 

        # making sure that d.official does not match any item in exogenous
        # the 'p' function gets the correct 'property' of that line
        eidlens = dict() 
        bw      = 0.0
        for d in self.fwdata :

            eidlens.update({ d : eidLen( d.entrez, d.organism )})
            if notNone([ r.match(d.official) for r in exo ]) :
                pass 
            elif d is self.bait : 
                pass 
            else : 
                bw += p(d)/eidlens[d]; 

        self.pseudo = 1.0 / bw / PSEUDO_LENGTH 

        for d in self.fwdata : 
                        
            counts_this_gene = 0.0 

            # dealing with background -- either single background count is exceeded
            # or peptide count decreased
            if self.background and self.background.get(d.official) : 

                background_exceeded=False 
                previous_fc=1000.0

                for fc in d.fxncounts : 

                    if background_exceeded or fc > self.background.get(d.official) \
                       or fc > previous_fc :
                        counts_this_gene += fc 
                        background_exceeded =   True
                    else : 
                        previous_fc = fc

            else : 
                counts_this_gene = p(d) 

            uselen = eidlens.get(d,PSEUDO_LENGTH) 
            d.setScore(counts_this_gene / bw / uselen) 


    def syncToEntrez( self, tryhard = True, debug = False, bestpepdb = 'RPHs', reference = hsgREF ) :

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if bestpepdb == 'RPMm':
            reference = mmgREF
        w  = 0
        W  = len(self.fwdata)
        for d in self.fwdata : 
            ( entrez, sym, org ) = desc_interpreter( d.desc, tryhard = tryhard, debug = debug,
                                                     bestpepdb = bestpepdb, reference = reference )

            d.setOfficial( sym )
            d.setEntrez( entrez )    
            d.setOrganism( org )


    def save(self,fname="",bait=None,concatenate=True,m2h=False,h2m=False,debug=False,no_zeros=True) :
        #nb that asMouse trumps asHuman

        def orgconvert(eid,converter) : 
            poseid=converter.get(eid,None) 
            if not poseid : 
                sys.stderr.write('Failed conversion for {}.\n'.format(eid)) 
                return
            elif type(poseid) is list : 
                weid=poseid[0] 
            else : 
                weid=poseid 

            if type(hsgREF[weid]) is list : 
                posref=hsgREF[weid][0] 
            else : 
                posref=hsgREF[weid] 

            worg=posref['Taxon'] 
            woff=posref['Symbol']

            return (weid,worg,woff) 

            
        if bait is None : 
            bait=self.bait 
            
        if not m2h and not h2m : 
            baitSym=bait.official 
            baitOrg=bait.organism 
            baitEnt=bait.entrez
        elif h2m : 
            try : 
                (baitEnt,baitOrg,baitSym)=orgconvert(bait.entrez,hs2mm) 
            except (TypeError,KeyError) : 
                baitSym=bait.official 
                baitOrg=bait.organism 
                baitEnt=bait.entrez
        elif m2h : 
            try : 
                (baitEnt,baitOrg,baitSym)=orgconvert(bait.entrez,mm2hs) 
            except (TypeError,KeyError) : 
                baitSym=bait.official 
                baitOrg=bait.organism 
                baitEnt=bait.entrez

        if concatenate : 
            sys.stderr.write('WARNING: if you scored by a non-additive measure, this will wreck things.\n') 
            outData=list() 
            addedData=set() 
            dataMid=dict() 

            for d in self.fwdata : 

                if d.official not in addedData : 
                    dataMid.update({ d.official : d })    
                    addedData.add(d.official) 
                else : 

                    sys.stdout.write('Concatenating observations {} and {} (symbol {}).\n'.\
                     format(d.idn,dataMid[d.official].idn,d.official)) 
                   #for i in range(0,len(d.fxncounts)) : 
                   #    dataMid[d.official].fxncounts[i] += d.fxncounts[i] 

                   #dataMid[d.official].totalcounts=sum(dataMid[d.official].fxncounts) 
                   #dataMid[d.official].maxcount=max(dataMid[d.official].fxncounts) 

                    dataMid[d.official].score += d.score
                    dataMid[d.official].totalcounts += d.totalcounts 

            outData=list(dataMid.values()) 

        else : 
            outData=self.fwdata 



        if ( not fname ) :
            outfile=sys.stdout
        else :
            outfile=open(fname,"w") 

        outfile.write("ID\tBAIT\tPREY\tNSAF\tORGANISM-BAIT\tORGANISM-PREY\tENTREZ-BAIT\tENTREZ-PREY\tDATASET\tNOTES\n") 

        for d in outData :

            if no_zeros and d.score == 0.0 : 
                sys.stdout.write('  Skipping {} (score of 0)\n'.format(d.official)) 
                continue 

            if not m2h and not h2m : 
                woff=d.official 
                worg=d.organism 
                weid=d.entrez 
            elif h2m : 
                # we are converting genes from whatever TO MOUSE
                woff='' 
                worg='' 
                weid='' 

                try : 
                    (weid,worg,woff)=orgconvert(d.entrez,hs2mm) 
                except (TypeError,KeyError) : 
                    woff=d.official 
                    worg=d.organism 
                    weid=d.entrez 

            elif m2h : 
                # we are converting genes from whatever TO HUMAN
                woff='' 
                worg='' 
                weid='' 

                try : 
                    (weid,worg,woff)=orgconvert(d.entrez,mm2hs) 
                except (TypeError,KeyError) : 
                    woff=d.official 
                    worg=d.organism 
                    weid=d.entrez 

            notestring = 'prey:'+woff+'_len_'+repr(eidLen( weid, worg ))+'_raw_'+repr(d.totalcounts) 

            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
             self.name+'_'+str(d.idn),baitSym,woff,d.score,baitOrg,worg,baitEnt,weid,self.name,\
             notestring)); 

        outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
         self.name+'_PSEUDO',baitSym,'PSEUDO',self.pseudo,baitOrg,'00',baitEnt,'00',self.name,\
         'pseudocount')); 

        outfile.close() 

    def setBait(self,baitOfficial) :

        for d in self.fwdata : 
            if d.official==baitOfficial :
                self.bait=d  
                break 
        else : 
            sys.stderr.write('Bait could not be identified!') 
            self.bait=None 
    


#wildly inefficent but occasionally necessary
def seqRetter(descline) :

    tremblre=re.compile(r'.*tr\|([^\|]*)\|.*') 
    swisspre=re.compile(r'.*sp\|([^\|]*)\|.*') 

    if swisspre.match(descline) : 
        spid=swisspre.match(descline).group(1)
    elif tremblre.match(descline) : 
        spid=tremblre.match(descline).group(1)

    if not spid : 
        sys.stderr.write("\nMSpreprocess.SeqRetter: could not parse accession from {}.\n".format(descline))
        return

    call(['seqret','swissprot:'+spid,'.tmp.txt'])

    tmp=file(".tmp.txt","r")
    srdesc=tmp.readline()
    tmp.close()

    return srdesc 

def desc_interpreter( desc, tryhard = True, debug = False, bestpepdb = 'RPHs', reference = hsgREF ) : 

    def refertouser( faildesc ) :

        if faildesc in rb.dup : 
            if debug : 
                sys.stdout.write('Description in rbase.dup.\n')
            return rb.dup[ faildesc ] 

        sys.stderr.write("\nDATASET: Provided description \n\t'{}'\n  cannot be mapped to any entrez entry.".\
         format(faildesc ));

        while True : 
            sys.stderr.write("Please input the correct symbol manually> ")
            sys.stderr.flush() 
            fetchedSym=sys.stdin.readline().strip()
            if mu.yesorno("Use {} as a gene symbol?".format(fetchedSym)) :
                break 

        try :
            if type(reference['symbol'][fetchedSym]) is list : 
                (e,i,o)  = referwithlist(faildesc,reference[fetchedSym]) 
                sys.stderr.write('Assigning entry to {}:{} (taxon {})\n'.format(e,i,o))
                sys.stderr.flush()
                return ( e, i, o )
            else : 
                sys.stderr.write('Assigning entry to {}:{} (taxon {})\n'.format(reference['symbol'][fetchedSym]['eid'],fetchedSym,\
                 reference['symbol'][fetchedSym]['taxid']))
                sys.stderr.flush()
                return( reference['symbol'][fetchedSym]['eid'], fetchedSym, reference['symbol'][fetchedSym]['taxid'] )  
        except (KeyError,TypeError) : 
            sys.stderr.write('Symbol {} doesn\'t match any records. Assigning EID and taxa 00.\n'\
            .format(fetchedSym))
            sys.stderr.flush()
            return( '00', fetchedSym, '00' )  

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def referwithlist( desc, possreclist ) :
        global rbase 
        sys.stderr.write("\nDATASET: Provided description\n\t{}\n is ambiguous.\n".format(desc))
        sys.stderr.write("Please choose an entry manually from the menu below or -1 to type in another symbol.\n");

        # possreclist is now actual RECORDS (type dicts) rather than symbols

        for w in range(0,len(possreclist)) : 
            sys.stderr.write("{: >4}----{:->16}----{:->8}----{:->8}\n".format(\
             w,possreclist[w]['eid'],possreclist[w]['symbol'],possreclist[w]['taxid']))

        sys.stderr.write("> ")
        userin = sys.stdin.readline()
        while True :
            k  = 0 
            try :
                k = int(userin.strip())
                if ( k == -1 ) :
                    return refertouser(desc)
                else:
                    sym    = possreclist[k]['symbol']
                    entrez = possreclist[k]['eid']
                    org    = possreclist[k]['taxid']
                return (entrez,sym,org)
            except (IndexError,ValueError) :
                sys.stderr.write("Invalid integer.\n> ")
                userin     = sys.stdin.readline()
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # now that entrez ids are a thing, we need to know first whether we're looking at 
    # an NCBI or a uniprot entry

    sym     = None 
    entrez  = None 
    org     = None 

    if debug:
        print( 'description to match: ' + desc + '\n' )
    
    me      = entrezre.match( desc )
    ms      = swisspre.match( desc ) 
    mt      = tremblre.match( desc ) 
    mg      = symbolre.match( desc ) 
    mp      = proteinre.match( desc )  

    froment = list() 
    fromsym = list() 
    fromswp = list() 
    fromtre = list() 
    fromsyn = list() 
    frompep = list() 

    if me : 
        froment = reference['eid'].get( me.group(1), '' ) 
    if mg : 
        fromsym = reference['symbol'].get( mg.group(1), '' ) 
    if ms : 
        fromswp = reference['swissprot'].get(ms.group(1).upper(),'') 
    if mt : 
        fromtre = reference['trembl'].get(mt.group(1).upper(),'')
    if mg : 
        fromsyn = reference['synonym'].get(mg.group(1),'') 
    if mp : 
        frompep = reference['peptide'].get(mp.group(1),'')

    mid     = dict() 
    fromall = list() 

    # concatenate all possible entries found
    for caught in [ froment, fromsym, fromsym, fromswp, fromtre, fromsyn, frompep ]  :
        # fromsym appears twice deliberately -- the symbol should be given
        # greater credence than synonyms
        if caught and type(caught) is list : 
            fromall.extend(caught) 
        elif caught and type(caught) is dict : 
            fromall.append(caught) 

    if debug : 
        sys.stderr.write("DEBUG:  Length of all matched regexes is {}.\n'".format(len(fromall)))

    for entry in fromall :
        if entry and entry['eid'] not in mid : 
            mid.update({ entry['eid'] : 1 }) 
        elif entry :
            mid[entry['eid']] += 1; 

    if mid : 
        choiceEnt = list() 
        choices   = list() 
        best      = max( mid.values() ) 

        for entryeid in mid : 
            if mid[entryeid] == best :
                choiceEnt.append(entryeid) 

        for ce in choiceEnt : 

            if type(reference['eid'][ce]) is list : 
                choices.extend(reference['eid'][ce]) 
            else : 
                choices.append(reference['eid'][ce]) 

        if debug : 
            sys.stderr.write('\nDEBUG: AMBIGUITY :\n') 

            sources = [('entrez', froment), ('symbol',fromsym), ('swiss',fromswp),
                       ('trembl', fromtre), ('synonym',fromsyn), ('peptide',frompep)] 

            for src,rec in sources : 
                if type(rec) is list : 
                    for r in rec : 
                        sys.stderr.write('    {:-<12}--{:->12}--{:->12}\n'.format(src,r['symbol'],r['eid'])) 
                elif rec:
                    sys.stderr.write('    {:-<12}--{:->12}--{:->12}\n'.format(src,rec['symbol'],rec['eid'])) 
            

            sys.stderr.write('    SCORING:\n') 
            for c in choices : 
                sys.stderr.write('    {:-<12}--{:-<12}..{:.>}\n'.format(c['symbol'],c['eid'],mid[c['eid']])) 

        if choices and type(choices) is list and len(choices) > 1 :
            (entrez,sym,org) = referwithlist(desc,choices); 

        elif choices and type(choices) is list :
            sym    = choices[0]['symbol'] 
            entrez = choices[0]['eid'] 
            org    = choices[0]['taxid'] 
        elif tryhard and ms : 
            try : 
                swacc,seq = E.fetchSw(ms.group(1),asTuple=True) 
                pepacc    = pep.bestPep(seq,db=bestpepdb) 
                sym       = reference['peptide'][pepacc]['Symbol'] 
                entrez    = reference['peptide'][pepacc]['EID'] 
                org       = reference['peptide'][pepacc]['Taxon'] 
            except (KeyError,ValueError,IOError) : 
                pass 
        elif tryhard and mt :
            try : 
                swacc,seq = E.fetchSw(mt.group(1),asTuple=True) 
                pepacc    = pep.bestPep(seq,db=bestpepdb) 
                sym       = reference['peptide'][pepacc]['Symbol'] 
                entrez    = reference['peptide'][pepacc]['EID'] 
                org       = reference['peptide'][pepacc]['Taxon'] 
            except (KeyError,ValueError,IOError) : 
                pass 
        elif tryhard and mp : 
            try :
                rec       = E.fetchPR(mp.group(1))
                sym       = E.PRgetSym(rec) 
                entrez    = E.PRgetEID(rec) 
                org       = E.PRgetOrg(rec) 
            except (KeyError,ValueError,IOError) : 
                pass 

        else:
            (entrez,sym,org) = refertouser(desc)
#            if mt and mt.group(1) and entrez != '00' and org != '00'  : 
#                rb.dup.update({ desc : (entrez,sym,org) }) 

    elif tryhard and ms : 
        try : 
            swacc,seq     = E.fetchSw(ms.group(1),asTuple=True) 
            pepacc        = pep.bestPep(seq,db=bestpepdb) 
            sym           = reference['peptide'][pepacc]['Symbol'] 
            entrez        = reference['peptide'][pepacc]['EID'] 
            org           = reference['peptide'][pepacc]['Taxon'] 
        except (KeyError,ValueError,IOError) :
            (entrez,sym,org) = refertouser(desc)
            rb.dup.update({ desc : (entrez,sym,org) }) 
    elif tryhard and mt : 
        try : 
            swacc,seq     = E.fetchSw(mt.group(1),asTuple=True) 
            pepacc        = pep.bestPep(seq,db=bestpepdb) 
            sym           = reference['peptide'][pepacc]['Symbol'] 
            entrez        = reference['peptide'][pepacc]['EID'] 
            org           = reference['peptide'][pepacc]['Taxon'] 
        except (KeyError,ValueError,IOError) :
            (entrez,sym,org) = refertouser(desc)
            rb.dup.update({ desc : (entrez,sym,org) }) 
    elif tryhard and mp : 
        try :
            rec=E.fetchPR(mp.group(1))
            if E.PRgetEID(rec) not in reference['eid'].keys() : 
                if debug :
                    sys.stderr.write("DEBUG:   remotely matched description :\n    {}\n"\
                     .format(desc) +\
                    "        but it could not be mapped to current reference.\n")
                pass 
            else : 
                sym    = E.PRgetSym(rec) 
                entrez = E.PRgetEID(rec) 
                org    = E.PRgetOrg(rec) 

                if debug :
                    sys.stderr.write("DEBUG:   remotely matched description :\n    {}\n"\
                     .format(desc) +\
                    "        to EID : {} Symbol : {}\n".format(entrez,sym))

        except (KeyError,ValueError,IOError) : 
            pass 
    else :
        (entrez,sym,org) = refertouser(desc)
        #rbase.dup.update({ desc : (entrez,sym,org) }) 
        
    if not sym or not entrez or not org :
        (entrez,sym,org) = refertouser(desc)

    # an initial if/return shoud have caught cases where 'desc' is already in rbase.dup
    if entrez == '00' or entrez in reference['eid'].keys() :
        pass
        #        rbase.dup.update({ desc : (entrez,sym,org) }) 

    return (entrez,sym,org)

def save_dup() :
    rb.update_dup( rb.dup )


def eidLen( eid, org, suppress = True ) : 

    org = int(org)
    eid = int(eid) # some eids are '00' !
    
    if eid > 0 and eid in gREF[ org ]['eid']:

        possPeps = gREF[ org ][ 'eid' ][ eid ]['peptide'].split(";")
        isnp     = lambda acc : 'NP_' in acc and acc in pREF[org]
        isxp     = lambda acc : 'XP_' in acc and acc in pREF[org]

        if   any([ isnp(x) for x in possPeps ]) :
            return float("{0:.1f}".format(np.mean([ pREF[org][x]['len'] for x in possPeps if isnp(x) ])))
        elif any([ isxp(x) for x in possPeps ]) :
            return float("{0:.1f}".format(np.mean([ pREF[org][x]['len'] for x in possPeps if isxp(x) ])))


    if not suppress : 
        sys.stderr.write('NOTE: No satisfactory peptide accessions for eid {}, using average length of 375.\n'.format(eid)) 

    return PSEUDO_LENGTH

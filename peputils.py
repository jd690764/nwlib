import sys
from subprocess import check_output as co
from subprocess import CalledProcessError 
from os import environ 

CDTRACK=environ['BLASTDB'] ;
def tocsl(acclist) : 
    out='' ; 
    for acc in acclist[1:len(acclist)] : 
        if acc is not None : 
            out += acc + ',' ; 

    return out ; 
        


def pepfetch(accession,db='RPHs') : 

    try : 
        if type(accession) in [list,tuple,set] : 
            out=co(['blastdbcmd','-db',db,'-dbtype','prot','-entry',tocsl(accession),'-outfmt',\
            "%s"]).split(); 
        else : 
            out=co(['blastdbcmd','-db',db,'-dbtype','prot','-entry',accession,'-outfmt',\
            "%s"]); 
        return out ; 

    except CalledProcessError : 
        return None ; 


def peplen(accession,db='RPHs') : 

    try : 
        if accession is None :
            return None ; 

        if type(accession) in [list,tuple,set] : 
            out=list() ;
            for acc in accession : 
                out.append(int(co(['blastdbcmd','-db',db,'-dbtype','prot','-entry',tocsl(accession),'-outfmt',\
                "%l"]).strip())) ; 
        else:
            out=int(co(['blastdbcmd','-db',db,'-dbtype','prot','-entry',accession,'-outfmt',\
            "%l"]).strip()) ; 

        return out ; 
    except CalledProcessError : 
        return None ; 


def peptitle(accession,db='RPHs') : 
    try : 
        if type(accession) in [list,tuple,set] : 
            out=co(['blastdbcmd','-db',db,'-dbtype','prot','-entry',accession,'-outfmt',\
            "%t"]).strip().split('\n') ; 
        else : 
            out=co(['blastdbcmd','-db',db,'-dbtype','prot','-entry',accession,'-outfmt',\
            "%t"]).strip() ; 
        return out ; 
    
    except CalledProcessError : 
        return None ; 


def pepfasta(accession,db='RPHs',asTuple=False) : 
    # this one is NOT ok with lists

    try:
        if asTuple : 
            out=list(co(['blastdbcmd','-db',db,'-dbtype','prot','-entry',accession,'-outfmt',\
            "%f"]).strip().split('\n')) ; 

            defline=out[0] ; 
            seq=""

            for line in out[1:len(out)] : 
                seq += line.strip() ; 

            return (defline,seq) ; 

        else : 
            out=co(['blastdbcmd','-db',db,'-dbtype','prot','-entry',accession,'-outfmt',\
            "%f"]) ; 

            return out ; 
    except CalledProcessError : 
        return None ; 

def bestPep(seq,db='RPHs') :

    return co('echo "' + seq + '" | blastp -db ' + db + ' -evalue 0.01 -outfmt "6 sacc" -query /dev/stdin | head -n1',shell=True).decode('utf-8').strip() ; 
    #return co(['echo','"',seq,'"','|','blastp','-db',db,'-outfmt','"6 sacc"',\
         #'-num_alignments','1','-query','/dev/stdin'],shell=True) ; 

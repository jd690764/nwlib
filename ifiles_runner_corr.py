import sys
import os.path
import pickle
import numpy as np
from lib.markutils import tabulate
from scipy.stats.distributions import norm
from statsmodels.sandbox.stats.multicomp import multipletests

IFILES = '/mnt/msrepo/ifiles/'

def filesForBkgd( infile ):

    alldatadict = dict() ;
    with open(sys.argv[1]) as f : 
        for line in f : 
            alldatadict.update({line.strip():dict()})
            
    return alldatadict

def mad(series) : 
    return np.percentile(np.abs(series-np.percentile(series,50)),50) ; 

def filterExptsByPseudoCountDistr( ddict ):
    
    # remove experiments where the pseudocount is high
    # relative to the other pseudocounts
    pseudodict     = { k : ddict[k]['PSEUDO'] for k in ddict }
    pskeys         = list(pseudodict.keys())
    pslogvals      = np.log10(list(pseudodict.values()))
    pslogmad       = mad(pslogvals) ; 
    pslogmedian    = np.percentile(pslogvals,50)
    pslvps_hi      = 1-norm.cdf((pslogvals-pslogmedian)/pslogmad)
    rejected_ds_hi = multipletests( pslvps_hi, alpha=0.05 )[0]

    # return data in a dictionary
    filteredExpts  = {  pskeys[i] : rejected_ds_hi[i] for i in range(len(pskeys))}
    return filteredExpts
    

def createBkgdFile( alldatadict, outfile ):

    
    for fn in alldatadict.keys() :
        f = open( IFILES + fn) ; 
        f.readline() ;
        for ll in tabulate(f) : 
            alldatadict[fn].update({ ll[2] : np.float(ll[3]) })
        f.close() ;

    filteredExpts  = filterExptsByPseudoCountDistr( alldatadict )        
    allfns         = list() ;
    allsyms        = set() ;
    
    for expt in filteredExpts : 
        if filteredExpts[ expt ]: 
            sys.stderr.write('WARNING: '+ expt +' is a pseudocount outlier '+\
                             'and as such will NOT be incorporated\n')
        else : 
            allfns.append(expt) ; 
            allsyms |= set(alldatadict[expt].keys())

    allsyms        = list(allsyms) # all gene symbols in bkgd file 
    grid           = np.zeros((len(allsyms),len(allfns)))

    # fill up the matrix 
    for i in range(len(allsyms)) :  
        for j in range(len(allfns)) : 
            grid[i][j] = alldatadict[allfns[j]].get(allsyms[i],alldatadict[allfns[j]]['PSEUDO'])

    loggrid    = np.log10(grid) ;
    logmads    = np.fromiter( ( mad(loggrid[i,:]) for i in range(len(allsyms)) ),dtype=np.float) ; 
    logmeans   = np.fromiter( ( np.mean(loggrid[i,:]) for i in range(len(allsyms)) ),dtype=np.float) ; 
    logmedians = np.fromiter( ( np.percentile(loggrid[i,:],50) for i in range(len(allsyms)) ),dtype=np.float) ; 

    with open( outfile,'wb')as f : 
        pickle.dump( (allsyms,allfns,loggrid), f )


def mkBkgdFile( infile, outfile ):
    files = filesForBkgd( infile )
    createBkgdFile( files, outfile )

if __name__ == "__main__":

    if len(sys.argv[:]) < 3 : 
        exit ; 

    mkBkgdFile( sys.argv[1], sys.argv[2] )
    

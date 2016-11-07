import os.path
from lib import config as cf
import re

def makeGct( ifiles, outfile, scoreCol = 3 ):
    # assemble a gct file from a list of ifiles using
    # preys as genes, the filenames as the experiment and
    # the nsaf as the value (unless score says column 9
    # - then use counts

    expts  = list()
    genes  = dict()
    values = dict()
    
    for f in ifiles:
        if os.path.isfile(cf.ifilesPath + f):
            with open( cf.ifilesPath + f, 'rt' ) as fh:
                expt = re.match( '^(.+)_?\d{0,}\.i$', f ).group(1)
                expts.append( expt )
                for line in fh:
                    l     = line.strip().split('\t')
                    score = l[3]
                    if scoreCol == 9:
                        l[9] = re.sub(r'^.+_(\d+)$', r'\1', l[9])
                        score = l[9]

                    if expt not in values:
                        values[ expt ] = dict()
                    if l[2] not in genes:
                        genes[ l[2] ] = 1

                    if l[2] not in genes:
                        genes[ l[2] ] = 1

                    values[expt][l[2]] = score


    with open( outfile, 'wt' ) as outh:
        outh.write( '#1.2' + '\n' )
        outh.write( str(len(genes)) + '\t' + str(len(expts)) + '\n' )
        outh.write( 'gene' + '\t' + '\t'.join( expts ) + '\n' ) 
        for gene in genes:
            vs = list()
            vs.append( gene )
            for expt in expts:
                if gene in values[ expt ]:
                    v = values[ expt ][gene]
                else:
                    v = ''
                vs.append( v )
            outh.write('\t'.join( map( str, vs )) + '\n' )

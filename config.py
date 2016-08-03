
djPath          = '/usr/local/share/py/djscripts/'

superfamilyFtp  = 'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links'
cddIdTableFtp   = 'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz'

uniprotUrl      = 'http://www.uniprot.org/uniprot/'


controlFiles    = '/mnt/reference/control_panels/'

ifilesPath      = '/mnt/msrepo/ifiles/' 

properNounsPath = '/home/mrkelly/Lab/errata/propers.txt'

referencePath   = '/usr/local/share/py/djscripts/data/pickles/'

publicDataPath  = '/mnt/reference/'
biogridPath     = publicDataPath + 'biogrid_latest'
#emiliomePath    = publicDataPath + 'complexes/emiliome.i'
emiliomePath    = 'data/pickles/emiliome.i'
emiliomeV2Path    = 'data/pickles/emiliomeV2.i'
bioplexPath     = publicDataPath + 'bioplex.i'
preppiPath      = publicDataPath + 'preppi_150727_lr600.i'
domainPairsPath = publicDataPath + 'bioplex_predicted.tsv'


filesDict       = { 'hsg' : referencePath + 'hsg_latest' ,
                    'mmg' : referencePath + 'mmg_latest' ,
                    'hmg' : referencePath + 'hsmmg_latest' ,
                    'hsp' : referencePath + 'hsp_latest' ,
                    'mmp' : referencePath + 'mmp_latest' ,
                    'cdd' : referencePath + 'cdd_latest' ,
                    'h2m' : referencePath + 'h2m_latest' ,
                    'm2h' : referencePath + 'm2h_latest' ,
                    'dup' : referencePath + 'dup_latest' ,
}

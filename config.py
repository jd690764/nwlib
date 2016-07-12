

superfamilyFtp  = 'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links'
cddIdTableFtp   = 'ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz'

uniprotUrl      = 'http://www.uniprot.org/uniprot/'


controlFiles    = '/mnt/reference/control_panels/'

ifilesPath      = '/mnt/msrepo/ifiles/' 

properNounsPath = '/home/mrkelly/Lab/errata/propers.txt'

referencePath   = '/usr/local/share/py/djscripts/data/pickles/'

publicDatapath  = '/mnt/reference/'
biogridPath     = publicDatapath + 'biogrid_latest'
emiliomePath    = publicDatapath + 'complexes/emiliome.i'
bioplexPath     = publicDatapath + 'bioplex.i'
preppiPath      = publicDatapath + 'preppi_150727_lr600.i'
domainPairsPath = publicDatapath + 'bioplex_predicted.tsv'


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

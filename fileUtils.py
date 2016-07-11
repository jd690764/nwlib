import zipfile
import os
import pycurl, sys, traceback, logging
import gzip
import yaml

def unzip(zipFilePath, destDir):
    zfile = zipfile.ZipFile(zipFilePath)
    for name in zfile.namelist():
        (dirName, fileName) = os.path.split(name)
        if fileName == '':
            # directory
            newDir = destDir + '/' + dirName
            if not os.path.exists(newDir):
                os.mkdir(newDir)
        else:
            # file
            fd = open(destDir + '/' + name, 'wb')
            fd.write(zfile.read(name))
            fd.close()
    zfile.close()

def downloadFromUrl( url, fpath ):

    with open( fpath, 'wb') as fh:
        curl = pycurl.Curl()
        curl.setopt(pycurl.URL, url )
        curl.setopt(pycurl.FOLLOWLOCATION, 1)
        curl.setopt(pycurl.MAXREDIRS, 5)
        curl.setopt(pycurl.CONNECTTIMEOUT, 30)
        curl.setopt(pycurl.TIMEOUT, 300)
        curl.setopt(pycurl.NOSIGNAL, 1)
        curl.setopt(pycurl.WRITEDATA, fh )
        try:
            curl.perform()
        except:
            traceback.print_exc(file=sys.stderr)
            sys.stderr.flush()
            logging.debug('Wrote {} '.format(sys.argv[1]))
        curl.close()

def gunzip( gzFile, path, destFile ):
    if os.path.exists( gzFile ) and os.path.exists( path ):
        with gzip.open( gzFile, 'rt' ) as src:
            with open( destFile, 'wt' ) as dest:
                dest.writelines( src )
    else:
        print( 'check if paths exist!' )


def read_config( yamlfile ):

    rescue_fn    = ''
    outfilename  = ''
    youtfilename = ''

    # read in config data from yaml file
    infile       = open( yamlfile )
    yout         = yaml.load( infile.read() )
    infile.close() 

    ds_dicts     = yout['datasets']
    # these should have the fields : infilename qualify bait convert

    for dsd in ds_dicts : 
        lost_files = 0 
        if not os.path.isfile('./'+dsd['infilename']) and not os.path.isfile('/mnt/msrepo/ifiles/'+dsd['infilename']) : 
            print('File '+dsd['infilename']+' not found.') ;
            lost_files += 1 ; 

    if lost_files > 0 : 
        raise IOError

    ALPHA_HI     = yout['options'].get('alpha_hi',0.01)
    ALPHA_LO     = yout['options'].get('alpha_lo',0.05)
    rescue_deg   = yout['options'].get('degree_to_rescue',1) ;
    valid_quals  = yout['options'].get('valid_degree_quals',{'wt',}) ; 
    if type(valid_quals) is list : 
        valid_quals = set( valid_quals )

    nwd          = yout['options'].get('network-wide_degree',False) ;

    if yout.get('files') : 
        rescue_fn    = yout['files'].get('rescue','') ; 
        youtfilename = yout['files'].get('outfile','') ;
        yidbfilename = yout['files'].get('idb','') ; 

    public_dicts = yout['public'] ;
    # this should be a list of dicts
    # each dict (1/public file) should have the fields : infilename qualify convert misncore minweight
    # public datasets have NO baits, are NOT DIRECTED, and are NEVER compared to negative controls

    # get list of symbols to rescue
    if os.path.isfile(rescue_fn) : 
        rescue_f = open(rescue_fn,'r') ; 
    elif rescue_fn : 
        sys.stderr.write('Invalid path provided for rescue file.\n') ; 
        os._exit(1) ; 
    else : 
        rescue_f = None

    # figure out output file name
    if youtfilename :
        outfilename = youtfilename ; 
    else : 
        sys.stderr.write('No valid output file name provided!\n') ; 
        os._exit(1) ; 

    # sort out idb output file name
    if yidbfilename : 
        idbfilename = yidbfilename ;


    conf         = dict( )
    conf[ 'ds_dicts' ] = ds_dicts
    conf[ 'ALPHA_HI' ] = ALPHA_HI
    conf[ 'ALPHA_LO' ] = ALPHA_LO
    conf[ 'rescue_deg' ] = rescue_deg
    conf[ 'valid_quals' ] = valid_quals
    conf[ 'nwd' ] = nwd
    conf[ 'public_dicts' ] = public_dicts
    conf[ 'rescue_f' ] = rescue_f
    conf[ 'outfilename' ] = outfilename
    conf[ 'idbfilename' ] = idbfilename

    return conf

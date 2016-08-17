import re


REGEX_COMPILED_TYPE=type(re.compile("foo")) ;
# for type comparisons later on

VALID_FIELDS_INTERACTIONS=["interID", "system", "systemType", "Author", "pmid", "throughput", "score", "modification", "phenotypes", "qualifications", "tags", "srcDB"] 
VALID_FIELDS_INPUTS=["interID", "entrezA", "entrezB", "biogridA", "biogridB", "systematicA", "systematicB", "officialA", "officialB", "synonymsA", "synonymsB", "system", "systemType", "Author", "pmid", "organismA", "organismB", "throughput", "score", "modification", "phenotypes", "qualifications", "tags", "srcDB"] 
VALID_FIELDS_NODES=["entrez","biogrid","systematic","official","synonyms","organism"]

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
r"ACT[ABG].*",\
r"MY[OH].*",\
r'FLN[ABC].*',\
r'VIM.*',\
r'ruvbl.*',\
}

highly_variable={\
r'HIST.*',\
r'AKR.*',\
r'CSNK.*',\
}

pseudogenes={\
r'LOC.*',\
r'Gm[0-9]\+',\
}

exogenous_regexes=exores

ribosomal_regexes = { r"E[IET]F" } ;

# ribosomal (eukaryotic) initiation,elongation,and termination factors

specious_experiments = { "Affinity Capture-RNA",\
                         "Co-localization",\
                         "Dosage Growth Defect",\
                         "Dosage Lethality",\
                         "Dosage Rescue",\
                         "Negative Genetic",\
                         "PCA",\
                         "Phenotypic Enhancement",\
                         "Phenotypic Suppression",\
                         "Positive Genetic",\
                         "Protein-RNA",\
                         "Synthetic Growth Defect",\
                         "Synthetic Haploinsufficiency",\
                         "Synthetic Lethality",\
                         "Synthetic Rescue",\
                         "Two-hybrid"}

especially_rigorous_experiments = { "Biochemical Activity",\
                                    "Co-crystal Structure",\
                                    "Co-purification",\
                                    "Reconstituted Complex"}


class subfilter(object):

    def __init__(self, rxstr, attribute):

        if ( attribute not in VALID_FIELDS_INTERACTIONS ) and (attribute not in VALID_FIELDS_NODES ):
            sys.stderr.write("ERROR: {} is not a valid field of the interaction class.\n".format(attribute)) ;
            return ;
        else:
            self.rxstr=rxstr ;
            self.attribute=attribute ;
            self.regex=re.compile(rxstr,re.IGNORECASE) ;


    def test(self, query):

        if ( self.regex.match(query) ):
            return True ; 
        else:
            return False ;

class Filter(object):

    def __init__( self, includes = set(), excludes = set(), debug = False ) :
        self.includes = includes ;
        self.excludes = excludes ;
        self.debug    = debug ;


    def apply(self,queries):
        # multiple queries, multiple subfilters
        out=set() ;

        for q in queries :
            if ( self.test(q) ) :
                out.add(q) ;

        return out ;


    def test(self,query):
        # one query, multiple subfilters

        for i in self.includes :
            if ( not i.test(getattr(query,i.attribute)) ) :
                if (self.debug):
                    sys.stderr.write("DEBUG:Filter.test: {!r} fails inclusion match {}.\n".format(getattr(query,i.attribute),i.rxstr)) ;
                return False ;

        for e in self.excludes :
            if ( getattr(query,e.attribute) and e.test(getattr(query,e.attribute))) :
                if (self.debug) :
                    sys.stderr.write("DEBUG:Filter.test: {!r} satisfies exclusion match {}.\n".format(getattr(query,e.attribute),e.rxstr)) ;
                return False ;

        return True ;

def bg_regex_assembler_excl(iorg="") :
    return bg_regex_assembler(iorg)[0]

def bg_regex_assembler_incl(iorg="") :
    return bg_regex_assembler(iorg)[1]

def bg_regex_assembler(iorg="") :

    default_regexes=exores | endores | cytores | highly_variable | pseudogenes ; 

    exn=set() ;
    exi=set() ;

    for r in default_regexes :
        exn.add(subfilter(r,"official")) ;

    for r in ribosomal_regexes : 
        exn.add(subfilter(r,"official")) ;

    for r in specious_experiments :
        exi.add(subfilter(r,"system")) ;

    if ( iorg ):

        inn={subfilter(iorg,"organism")};

        return( Filter(excludes=exi),Filter(excludes=exn,includes=inn)) ;
    else : 

        return( Filter(excludes=exi),Filter(excludes=exn)) ;

def exogenous_regex_assembler(debug=False) : 
    from lib.rbase import exores as exogenous_regexes

    exn=set() ; 

    for r in exogenous_regexes : 
        exn.add(subfilter(r,"official")) ;

    return Filter(excludes=exn,debug=debug) ; 

def exogenous_and_isoformtastic(debug=False) :
    from lib.rbase import highly_variable as hv
    from lib.rbase import exores 
    exn=set() ; 

    for r in exogenous_regexes | hv : 
        exn.add(subfilter(r,"official")) ;

    return Filter(excludes=exn,debug=debug) ; 

def create_rasfilter() : 
    warn(DeprecationWarning('I (Mark) can\'t vouch for create_rasfilter right now 3/1/2016') )

    nsubfs=[ subfilter(r,'official') for r in rbase.exores ] ; 
    nsubfs.extend([ subfilter(r,'official') for r in rbase.endores ]) ; 

    #isubfs=[ subfilter(r'22939629','pmid') ] ; 
    isubfs=list() ;

    isubfs.extend([ subfilter(r,'system') for r in specious_experiments ]) ; 

    return (Filter(excludes=isubfs),Filter(excludes=nsubfs)) ; 



import tkinter,tkinter.filedialog
import sys
from lib import markutils as mu
import re
from os.path import isfile
from os import remove
import importlib
from django.forms import model_to_dict
from lib import rbase
from lib import filters


# note that field dictionaries still need to include things like entrezA. HOWEVER, some field information is exclusive to interactions and others are exclusive to nodes

# in different field times, which columns(in order) map to which fields of interaction
# or node data

fdms=( "interID" , "officialA" , "officialB" , "score" , "organismA" , \
       "organismB" , "entrezA" , "entrezB" , "srcDB", "tags" )


fd_biogrid=( "interID" , "entrezA" , "entrezB" , "biogridA" , "biogridB" , \
 "systematicA" , "systematicB" , "officialA" , "officialB" , "synonymsA" , \
 "synonymsB" , "system" , "systemType" , "Author" , "pmid" , "organismA" , \
 "organismB" , "throughput" , "score" , "modification" , "phenotypes" , \
 "qualifications" , "tags" , "srcDB" ) ; 

fd_emili=( "interID" , "entrezA" , "entrezB" , "officialA" , "officialB" ) ; 

fd_chem=("interID","biogridA","entrezA","systematicA","official","synonymsA","organismA","","",\
"Author","pmid","","biogridB","officialB","synonymsB",)
# resume, from biogrid chemical tsv

DEFAULT_FIELD_DICTIONARY =  fdms

mouse_taxid = '10090'
human_taxid = '9606'

FILTER      = { 'biogrid': filters.bg_regex_assembler,
                'bg_excl' : filters.bg_regex_assembler_excl,
                'bg_incl' : filters.bg_regex_assembler_incl,
                'exogeneous': filters.exogenous_regex_assembler,
                'isoform': filters.exogenous_and_isoformtastic }

keyer=lambda x : x.offical + '_' + x.entrez ;
# function to generate node 'key'
# almost always <symbol>_<entrez id>
# this key is used to look up nodes in node-containing dicts or for
# easy comparison between nodes even if they aren't mapped to the same
# memory object. For example, KRAS_3845 will refer to human KRAS
# in two different datasets even if they are separate dataset objects

FORCE_MATCH_QUAL=True
# I'm not sure why I would ever want this to be false tbh.

#debugout=open('.interactors_dbg_log.txt','w') ; 

def ee(ek1,ek2) : 
    """
        Using two edge keys e1 and e2, tests if they are equivalent.
        DIRECTED edges are equivalent only if their keys match.
        UNDIRECTED edges are equivalent if the contained node keys match.
    """

    if ek1 == ek2 : 
        return True ;
    else : 
        nk1_1 = ''
        nk2_1 = ''
        thequal='' ; 

        ek1i=iter(ek1)
        c=next(ek1i) ; 
        while c not in {'>','^'} : 
            nk1_1 += c ; 
            c=next(ek1i) ; 

        edgeicon = c; 
        if edgeicon == '>' : 
            return False ; 
            # means that this is a directed edge, so swapping
            # node key order shouldn't produce a match

        c=next(ek1i)
        # move past edge direction signifier
        while c != ':' : 
            nk2_1 += c ; 
            c=next(ek1i)

        try : 
            c=next(ek1i) ; 
            while True : 
                thequal += c ; 
                c=next(ek1i) ; 
        except StopIteration :
            pass ; 

        if nk2_1 + edgeicon + nk1_1 + ':' + thequal  == ek2 : 
            return True
        else : 
            return False ; 

def ei(ek1) :
    """
        Given an an edge key, return the edge key corresponding
        to an inverted edge. For undirected edges, this returns
        the input.
    """


    # inverts edge if directed, otherwise returns input
    nk1_1 = ''
    nk2_1 = ''
    thequal='' ; 

    ek1i=iter(ek1)
    c=next(ek1i) ; 
    while c not in {'>','^'} : 
        nk1_1 += c ; 
        c=next(ek1i) ; 

    edgeicon = c; 
    if edgeicon == '>' : 
        return ek1 ; 
        # means that this is a directed edge, so swapping
        # node key order shouldn't produce a match

    c=next(ek1i)
    # move past edge direction signifier
    while c != ':' : 
        nk2_1 += c ; 
        c=next(ek1i)

    try : 
        c=next(ek1i) ; 
        while True : 
            thequal += c ; 
            c=next(ek1i) ; 
    except StopIteration :
        pass ; 

    return nk2_1 + edgeicon + nk1_1 + ':' + thequal ; 

class interaction(object):
    """
    Contains data corresponding to observations of protein interactions.
    Fields are basically directly ripped from biogrid. Contains the following
    methods.

    interID---------a unique identifier for the observation
    author----------author of publication (biogrid only)
    directed--------whether the interaction is directed/not ...
    nodeA-----------pointer to node object. For directed nodes, this is the
                    "source." By my convention, in AP/MS the source node is
                    the bait protein.
    nodeB---------- pointer to the other node object.
    pmid------------pubmed ID of publication
    qualificiations-in non-biogrid data, used to mark different types of
                    interactions that should NOT be merged into the same
                    edge(e.g. interactions with mutant bait protein)
    score-----------in MOST cases, modified NSAF.
    srcDB-----------source database
    system----------experimental system
    systemType------physical or genetic (basically not used outside of biogrid)
    throughput------high-throughput or low-throughput (biogrid only)

    tags------------basically ignored
    modification----basically ignored
    phenotypes------basically ignored

    """
        
    def __init__(self,\
        # to A from B
        nodeA,
        nodeB,
        directed       = False,
        interID        = "",
        system         = "",
        systemType     = "",
        Author         = "",
        pmid           = "",
        throughput     = "",
        score          = "",
        modification   = "",
        phenotypes     = "",
        qualifications = "",
        tags           = "",
        srcDB          = "" ) :

        self.interID        = interID ; 
        self.Author         = Author ;
        self.directed       = directed ; 
        self.modification   = modification ;
        self.nodeA          = nodeA ;
        self.nodeB          = nodeB ;
        self.phenotypes     = phenotypes ;
        self.pmid           = pmid ;
        self.qualifications = qualifications ;
        self.score          = score ;
        self.srcDB          = srcDB ;
        self.system         = system ;
        self.systemType     = systemType ;
        self.tags           = tags ;
        self.throughput     = throughput ;
        
    def clone(self) :
        """
            Makes a copy of the interaction.
        """
        I = interaction(
            interID        = self.interID,\
            system         = self.system,\
            systemType     = self.systemType,\
            Author         = self.Author,\
            pmid           = self.pmid,\
            throughput     = self.throughput,\
            score          = self.score,\
            modification   = self.modification,\
            phenotypes     = self.phenotypes,\
            qualifications = self.qualifications,\
            tags           = self.tags,\
            srcDB          = self.srcDB,\
            nodeA          = self.nodeA,\
            nodeB          = self.nodeB)

        return I ;

    def __str__( self ):
        return( '; '.join( [
            'interID'        +': '+self.interID,
            'system'         +': '+self.system,
            'systemType'     +': '+self.systemType,
            'Author'         +': '+self.Author,
            'pmid'           +': '+str(self.pmid),
            'throughput'     +': '+self.throughput,
            'score'          +': '+str(self.score),
            'modification'   +': '+self.modification,
            'phenotypes'     +': '+self.phenotypes,
            'qualifications' +': '+self.qualifications,
            'tags'           +': '+self.tags,
            'srcDB'          +': '+self.srcDB,
            'nodeA'          +': '+str(self.nodeA),
            'nodeB'          +': '+str(self.nodeB)] ))
   #def set_nodeA(self,nodeA):
   #    self.nodeA=nodeA ;

   #def set_nodeB(self,nodeB):
   #    self.nodeB=nodeB ;

    def edgekey(self) : 
        """
            Generates key for an edge to which this interaction belongs.
        """
        if self.directed : 
            edgeicon='>' ; 
        else :
            edgeicon='^' ; 
        return self.nodeA.key + edgeicon + self.nodeB.key + ':' + self.qualifications

    def __contains__(self,thing) : 

        if thing in { self.nodeA,self.nodeB,self.nodeA.key,self.nodeB.key } : 
            return True ; 
        else : 
            return False ; 

class bgedge(object):
    """
        Container class of interactions. attributes:
        directed     - boolean if edge is directed or not.
        interactions - set of contained interactions
        key          - edge key. Takes the form <nk1>[^>]<nk2>:qual,
                       where the caret or wedge indicates an undirected or directed
                       edge respectively.
        meanscore    - mean score of contained interactions
        p            - p value for the interaction
        qual         - qualification string of contituent edges
        source       - interaction source string if container edges are uniform, otherwise.
        to           - destination node for interactions (same conventions as interaction class).
        totalscore   - total score of contained interactions
        weight       - number of interactions in edge
        whence       - sourcenode for interactions (i.e. APMS bait)
    """
    def __init__( self, interaction = None, directed = None, qual = '' ): 

        self.directed       = None
        self.interactions   = set()
        self.key            = ''
        self.meanscore      = 0.0
        self.p              = 1.0
        self.qual           = qual
        self.source         = '' 
        self.to             = None
        self.totalscore     = 0.0
        self.weight         = 0
        self.whence         = None 

        if interaction : 
            self.add_interaction( interaction ) ; 

    def add_interaction( self, interaction ):
        """
            bgedge.add_interaction(interaction)

            Adds the supplied interaction to this edge's substituent interactions,
            and updates class attributes accordingly.
        """

        if interaction in self.interactions : 
            return ; 

        if ( not self.to and not self.whence ) :
            self.whence         = interaction.nodeA ; 
            self.to             = interaction.nodeB ; 
        elif not { self.to , self.whence } ^ { interaction.nodeA, interaction.nodeB } : 
            # symmetric difference is empty ==> no overlap in nodes on this edge w/ incoming interaction
            pass ; 
        else : 
            raise KeyError('Edge\'s nodes are not empty, but incoming interaction\'s nodes do not match.\n') ; 

        if self.directed is None : 
            self.directed   = interaction.directed ; 
        elif self.directed != interaction.directed : 
            raise ValueError('Edge is directed and incoming interaction is not ,or vice versa.\n') ; 

        if self.qual != interaction.qualifications : 
            raise KeyError('Edge\'s and incoming interaction do not have matching qual/qualificaitons.\n') ; 

        self.interactions.add( interaction ) ; 
        self.weight += 1 ;
        try :
            self.totalscore += float(interaction.score) if interaction.score is not None else 0.0 ; 
            self.meanscore   = self.totalscore / self.weight ;
            
        except ValueError :
            pass ;

        if not self.key : 
            if self.directed  : 
                edgeicon='>' 
            else : 
                edgeicon='^'
            self.key = self.whence.key + edgeicon + self.to.key + ':' + self.qual

        self.determine_origins() ;

    def remove_interaction(self,interaction):
        """
            bgedge.remove_interaction(interaction)

            removes substituent interactions from edge,
            and updates attributes accordingly
        """

        if interaction not in self.interactions : 
            return ; 
        self.interactions.remove(interaction) ; 
        self.weight -= 1 ;
        try :
            self.totalscore -= float(interaction.score) ;
            self.meanscore   = self.totalscore / self.weight ; 
        except ValueError :
            pass ;

        self.determine_origins() ;


    def nkeys(self) :
        """
            returns a tuple of keys to the nodes
            at either end of this edge
        """
        return ( self.whence.key, self.to.key ) ;

    def determine_origins(self) :
        """
            changes source attribute depending on substituent
            interaction source strings.
        """

        istr = "" ;

        for i in self.interactions :
            if not istr : 
                istr = i.srcDB ;
            elif istr and i.srcDB != istr :
                istr = 'mixed' ;
                break ; 

        self.source = istr ;

    def connects(self,nodea) :
        """
            If the argument supplied is one of the nodes
            connected by this edge, returns the other node.

            Otherwise, returns None.
        """

        if nodea is self.whence :
            return self.to
        elif nodea is self.to : 
            return self.whence
        else : 
            return None ; 

    def connects_key(self,nodeakey) : 
        """
            If the argument supplied is a key to one of the nodes
            connected by this edge, returns the other node's key.

            Otherwise, returns None.
        """

        if nodeakey == self.whence.key :
            return self.to.key
        elif nodeakey == self.to.key : 
            return self.whence.key
        else : 
            return None ; 

    def sharesends(self, other, directed = False) :
        """
           returns true if:
               - undirected: the other edge connects the same nodes
               - directed: the other edge connects the same nodes in same orientation
           returns false otherwise
        """
        
        if {other.whence.key,other.to.key} != {self.to.key,self.whence.key} : 
            return False ;
        elif directed and self.to != other.to :
            return False ;
        else : 
            return True ;

   #def sameNode(node1,node2):
   #    if ( node1.official == node2.official and node1.organism == node2.organism ):
   #        return True ;
   #    else: 
   #        return False ;

    def __contains__(self,thing) : 

        if thing in { self.to,self.whence,self.to.key,self.whence.key } : 
            return True ; 
        else : 
            return False ; 

class node(object):
    """
        Contains node (gene) information.
        node.biogrid      =   ??????
        node.debug        =   flag for debug messages
        node.entrez       =   entrez ID of node if it corresponds to a real gene, otherwise 00
        node.key          =   key (symbol_entrez id)
        node.official     =   "official" gene symbol as used by entrez
        node.organism     =   taxid by NCBI conventions, 9606 for humans, 10090 for mice
        node.synonyms     =   "other genes" 
        node.systematic   =   biogrid relic
        node.interactions =   set containing interactions that link to this node
        node.edges        =   set containing edges that link to this node
        node.partners     =   set containing nodes that are linked by edges to this
    """

    def __init__( self, official, organism, key, entrez = "", biogrid = "", systematic = "", synonyms = "", debug = False ) :

        self.biogrid      = biogrid ;
        self.debug        = debug ; 
        self.entrez       = entrez ; 
        self.key          = key ; 
        self.official     = official ; 
        self.organism     = organism ; #this is the NCBI taxonomy id
        self.synonyms     = synonyms ; 
        self.systematic   = systematic ; 

        self.interactions = {} ; # key: biogrid id # ; value: interaction
        self.partners     = {} ; # key : partner official name + organism; value : node 
        self.edges        = {} ; # key : partner official name ; value: bgedge

    def __str__( self ):
        return( str(self.key) )
    
    def nneighbors( self, within_edge_set = None ) :
        """
            node.degree(self,within_edge_set = None) :

            returns the number of nodes that link to this one.
            If within_edge_set is supplied, ONLY count nodes
            that are linked by one or more of the provided edges

        """
        # python doesn't pass by reference
        if not within_edge_set : 
            return len([ k for k in list(self.partners.keys()) ])
        else : 
           return len([ k for k in list(self.partners.keys()) \
            if self.binds(k,within_edge_set=within_edge_set) ]) ; 
            # of nodes (because the edges dict is keyed by node) that have
            # some edges in the acceptable edge set (viz. second line does not produce empty set)
    
    def degree( self, within_edge_set = None ) :
        """
            node.degree(self,within_edge_set = None) :

            returns the number of nodes that link to this one.
            If within_edge_set is supplied, ONLY count nodes
            that are linked by one or more of the provided edges

        """
        # python doesn't pass by reference
        if not within_edge_set : 
            return len([ k for k in list(self.edges.keys()) ])
        else : 
            if type(next(iter(within_edge_set))) is str : 
               #return sum([ 1 if k in within_edge_set else 0  for k in self.edges.keys() ]) ;
               return len([ k for k in list(self.edges.keys()) \
                if { e.key for e in self.edges[k] } & within_edge_set  ]) ; 
            else : 
                return sum([ 1 if v & within_edge_set else 0 for v in list(self.edges.values()) ])
            # of nodes (because the edges dict is keyed by node) that have
            # some edges in the acceptable edge set (viz. second line does not produce empty set)

    def binds( self, othernode, within_edge_set = None ) : 
        """
            node.binds( self, othernode, within_edge_set = None ) :

            returns a boolean that indicates whether a node 
            (othernode) binds to this one.
            If within_edge_set is supplied, it ONLY considers
            edges that are part of the provided set

        """

        if within_edge_set is not None and not within_edge_set : 
            return False
        # handles the case where within_edge_set is an empty set

        if type(othernode) is str  : 
            if othernode not in list(self.edges.keys()) : 
                return False ; 
            elif within_edge_set : 
                if type(next(iter(within_edge_set))) is bgedge : 
                    return len( {e for e in self.edges[othernode] } & within_edge_set ) > 0 ;
                else : 
                    return len( { e.key for e in self.edges[othernode] } & within_edge_set) > 0
            else : 
                return True ; 

        elif type(othernode) == type(self) : 
             
            if othernode.key not in list(self.edges.keys()) :
                return False ;
            elif within_edge_set : 
                if type(next(iter(within_edge_set))) is bgedge : 
                    return len( {e for e in self.edges[othernode.key] } & within_edge_set ) > 0 ;
                else : 
                    return len({ e.key for e in self.edges[othernode.key] } & within_edge_set) > 0 ; 
            else : 
                return True ; 
        elif type(othernode) in {list,tuple,set} : 
            return [ self.binds(subnode,within_edge_set=within_edge_set) for subnode in othernode ]
        else : 
            raise TypeError('Argument to binds must be of type str or of type interactors.node\n') ; 


class dataSet(object):

    def __init__( self, nodes = None, default_organism = "9606", i_filter = None, n_filter = None, debug = False,
                  superdebug = False, default_src = "", correction_dict = None, interactions = None, edges = None ):

        self.default_organism   = default_organism ; 
        self.default_src        = default_src
        self.i_filter           = i_filter ; # interaction filter?
        self.n_filter           = n_filter ; # node filter? 
        self.infilenames        = list()  ; # filenames from which the network is constructed
        self.keys               = None ; # keys of nodes in this network
        self.superdebug         = superdebug ;
        self.debug              = True if self.superdebug else debug 

        if self.n_filter is not None and self.n_filter in FILTER:
            self.n_filter = FILTER[ self.n_filter ]()

        if self.i_filter is not None and self.i_filter in FILTER:
            self.i_filter = FILTER[ self.i_filter ]()
            
        # this will again be a set of interactions, but not necessarily with complete data
        if interactions is None :
            self.the_data = dict() ;
        else : 
            self.the_data = interactions ; # allows for a "starting" set of interactions

        if edges is None :
            self.edges    = dict() ;
        else : 
            self.edges    = edges ; # allows for a "starting" set of edges

        if nodes is None :
            self.nodes    = {} ; # empty dictionary
        else : 
            self.nodes    = nodes ; # allows for a "starting" dictionary
            self.keys     = set( self.nodes.keys() ) ; 

        if ( not correction_dict ) :
            self.correction_dict = dict() ;
        else : 
            self.correction_dict = correction_dict ; 

        #if ( self.debug ) :
            #debugout.write("DEBUG:   Sizes at initialization:\n Interactions : {} ; Nodes : {} ; Corrections : {} \n"\
            #.format(len(self.the_data),len(self.nodes),len(self.correction_dict))) ; 


    def _add_node( self, node ):

        """
            adds a new node to the network
            returns True if node was added or already exists, otherwise False
        """
        if not node.key in self.nodes.keys():
            if self.n_filter and not (self.n_filter and self.n_filter.test( node )):
                # the node is filtered out
                return False 
            else:
                if self.debug : 
                    debugout.write( "DEBUG:   Creating node {}.\n".format( node.key ))
                self.nodes.update( { node.key : node } )
                
        return True


    def _add_interaction( self, iaction ):

        """
            adds a new interaction to the network. Depepends on 
            self._add_node to add the nodes before this can be called
            returns True is the interaction was created/existed
               and False otherwise
        """

        if self.i_filter and not ( self.i_filter and self.i_filter.test( iaction )):
            #if self.debug:
            #    debugout.write("DEBUG: Interaction {} excluded after filtering.\n".format( iaction.interID))
            return False
        else:
            if not iaction.interID in self.the_data.keys():
                self.the_data.update({ iaction.interID : iaction }) ;
            
                # if the interaction belongs to a previously uncharacterized edge
                if iaction.edgekey() not in self.edges and ei( iaction.edgekey()) not in self.edges : 

                    # create the edge
                    newedge = bgedge( interaction = iaction,
                                      directed    = iaction.directed,
                                      qual        = iaction.qualifications)

                    # add it to the dataSet's dict of edges
                    self.edges.update({ newedge.key : newedge }) ;
                
                    # connect partners of this edge(may be unecessary) ; 
                    newedge.whence.partners.update({ newedge.to.key : newedge.to }) ; 
                    newedge.to.partners.update({ newedge.whence.key : newedge.whence }) ; 
                
                    if newedge.whence.key not in newedge.to.edges : 
                        newedge.to.edges.update({ newedge.whence.key : { newedge } }) ; 
                    else : 
                        newedge.to.edges[newedge.whence.key].add(newedge) ; 
                    
                    if newedge.to.key not in newedge.whence.edges : 
                        newedge.whence.edges.update({ newedge.to.key : { newedge } }) ; 
                    else : 
                        newedge.whence.edges[newedge.to.key].add(newedge) ; 

                # if the edge was previously characterized, but backwards
                elif ei( iaction.edgekey()) in self.edges : 
                    self.edges[ ei( iaction.edgekey())].add_interaction( iaction )
                else :
                    self.edges[ iaction.edgekey()].add_interaction( iaction )

                iaction.nodeA.interactions.update({ iaction.interID : iaction })
                iaction.nodeB.interactions.update({ iaction.interID : iaction })

        return True


    def _load_data_from_db( self, infobj, logger ):

        from django.core.exceptions import FieldDoesNotExist
        from network.models import Interaction as obj
        from django.db.models import Q
        import operator
        from functools import reduce
        
        infobj      = [ i.upper() for i in infobj ]
        
        # assemble filters
        kwargs_f    = { 'organisma': self.default_organism, 'organismb': self.default_organism }
        kwargs_f[ 'srcdb__in' ] = infobj

        # if we have nodes already, we are looking for interactions between those nodes...
        if self.nodes :
            eids   = [ k.split('_')[1] for k in self.keys ]
            q_list = [Q(x) for x in [ ('entreza__in', eids), ('entrezb__in', eids )]]
            
        kwargs_e    = {}
        alldata     = ''
    
        logger.debug('had q_list')
        if self.i_filter[ 'clude' ] == 'exclude' :
            logger.debug('exclude')
            kwargs_e[ self.i_filter[ 'column' ] + '__in' ] = self.i_filter[ 'values' ]
            alldata  = obj.objects.filter( reduce(operator.or_, q_list), **kwargs_f ).exclude( **kwargs_e ).values( 'entreza', 'entrezb' ) 
        elif self.i_filter[ 'clude' ] == 'include' :
            logger.debug('include')
            kwargs_f[ self.i_filter[ 'column' ] + '__in' ] = self.i_filter[ 'values' ] 
            alldata  = obj.objects.filter( reduce(operator.or_, q_list), **kwargs_f ).values( 'entreza', 'entrezb' ) 

        logger.debug( len(alldata) )
        #alldata = [ model_to_dict( data ) for data in alldata ]
        logger.debug( len(alldata) )
        eids    = [ [r['entreza'], r['entrezb']] for r in alldata ]
        logger.debug( len(eids) )
        eids    = list(set([ eid for elem in eids for eid in elem ]))
        logger.debug( len(eids) )

        kwargs_f[ 'entreza__in' ] = eids
        kwargs_f[ 'entrezb__in' ] = eids

        # ... but we also need interactions between neighbors
        if self.i_filter[ 'clude' ] == 'exclude' :
            logger.debug('exclude')
            kwargs_e[ self.i_filter[ 'column' ] + '__in' ] = self.i_filter[ 'values' ] 
            alldata  = obj.objects.filter( **kwargs_f ).exclude( **kwargs_e )
        elif self.i_filter[ 'clude' ] == 'include' :
            logger.debug('include')
            kwargs_f[ self.i_filter[ 'column' ] + '__in' ] = self.i_filter[ 'values' ] 
            alldata  = obj.objects.filter( **kwargs_f )
        
        # we still need to deal with n_filter somehow
        logger.debug( len(alldata) )
        return alldata
        
            
    def load( self, infobj, fd = None, h2m = False, m2h = False, directed = False, qualify = '',
              force_qualify = False, force_score = None, logger = '' ) :
        
        self.infilenames.append( infobj ) ;

        alldata     = self._load_data_from_db( infobj, logger )
        added_iKeys = set(self.the_data.keys()) ;
        added_nKeys = set(self.nodes.keys()) ;

        if fd :
            fields_dict = fd ;
        else :
            fields_dict = DEFAULT_FIELD_DICTIONARY ; 
        i = 1
        for data in alldata:
            data   = model_to_dict( data )
            node_was_filtered = False ;
            #logger.debug( data )
            if force_score is not None and data.get('score') is not None : 
                data['score'] = force_score ; 

            #if not i % 100:
            #    logger.debug( i )
            i += 1
            i_keys   = {"a" : "" , "b" : "" }  ; 

            ## nodes must be created before interactions
            for c in ( 'a','b') :

                sym         = data.get( 'symbol' + c )
                entrez      = data.get( 'entrez' + c )
                i_keys[c]   = sym + "_" + str(entrez) ; 

                new_node    = node( entrez       = entrez,
                                    biogrid      = data.get("biogrid"+c),
                                    systematic   = data.get("symbol"+c),
                                    official     = sym,
                                    synonyms     = data.get("synonyms"+c),
                                    # establish use of default organisms where not provided
                                    organism     = data.get( 'organism' + c, self.default_organism ),
                                    key          = i_keys[c] ,
                                    debug        = self.superdebug) ;

                if h2m : 
                    new_node  = mouseify_node(new_node) ; 
                    i_keys[c] = new_node.key ; 

                elif m2h : 
                    new_node  = humanize_node(new_node)
                    i_keys[c] = new_node.key ;

                if ( new_node.key not in added_nKeys ):
                    node_was_filtered = not( self._add_node( new_node ))
                    if not node_was_filtered:
                        added_nKeys.add(new_node.key)
                    
            if ( node_was_filtered ) :
                continue ; 

            qualify_line = qualify[ data.get( 'srcdb' )]
            if force_qualify : 
                thequal = qualify_line
            else : 
                thequal = data.get( 'qualifications', qualify_line ) ; 
            #logger.debug( 'b: "' + i_keys['b']+'"' )
            #logger.debug( self.nodes )
            if data["interid"] not in added_iKeys :
                thisinteraction = interaction( interID        =    data.get("interid"),\
                                               system         =    data.get("system"),\
                                               systemType     =    data.get("systemtype"),\
                                               Author         =    data.get("author", ''),\
                                               pmid           =    data.get("pmid"),\
                                               throughput     =    data.get("throughput"),\
                                               score          =    data.get("score", ''),\
                                               modification   =    data.get("modification", ''),\
                                               phenotypes     =    data.get("phenotypes", ''),\
                                               qualifications =    thequal,\
                                               tags           =    data.get("tags", ''),\
                                               directed       =    directed,\
                                               srcDB          =    data.get("srcdb"),\
                                               nodeA          =    self.nodes[i_keys['a']],\
                                               nodeB          =    self.nodes[i_keys['b']],) ; 
            else:
                logger.debug("WARNING: Interaction {} on line {} has already been added!\n".format(data["interid"],i))
                continue ;

            if self._add_interaction( thisinteraction ):
                added_iKeys.add(data["interID"])
                
        logger.debug("Completed.") ;
        
    def parse( self, infobj, sep = '\t', fd = None, h2m = False, m2h = False, directed = False, qualify = '',
               force_qualify = False, force_score = None ) :

        if type(infobj) is str : 
            infobj=open(infobj) ;

        self.infilenames.append( infobj.name )
        
        filelines   = sum( 1 for line in infobj ) ;
        sys.stderr.write('Parsing {} ({} lines).\n'.format(infobj.name,filelines)) ; 
        infobj.seek(0) ; # return pointer to beginning of file

        added_iKeys = set(self.the_data.keys()) ;
        added_nKeys = set(self.nodes.keys()) ;
        headerline  = infobj.readline() ;
        headers     = headerline.split(sep) ;
        fields_dict = fd if fd else DEFAULT_FIELD_DICTIONARY
        dataline    = infobj.readline()
        i           = 1
        converted   = dict()
        
        while dataline : 

            # skip all commented lines ; first uncommented line is header
            if dataline[0] == '#' : 
                continue
            
            node_was_filtered = False
            dataline          = dataline.strip()
            vals              = dataline.split( sep )
            data              = dict(list(zip( fields_dict, vals)))
            i_keys            = { "A" : "" , "B" : "" }

            if force_score is not None and data.get('score') is not None : 
                data['score'] = force_score ; 

            ## nodes must be created before interactions
            for c in ( 'A','B') :

                sym         = data.get( 'official' + c )
                entrez      = data.get( 'entrez' + c )
                ori_key     = sym + "_" + entrez
                i_keys[c]   = ori_key

                new_node    = node(entrez       = entrez,
                                   biogrid      = data.get( "biogrid"+c ),
                                   systematic   = data.get( "systematic"+c ),
                                   official     = sym,
                                   synonyms     = data.get( "synonyms" + c ),
                                   organism     = data.get( 'organism' + c, self.default_organism),
                                   key          = ori_key,
                                   debug        = self.superdebug) ;
                if h2m :
                    if ori_key in converted:
                        new_node = converted[ ori_key ]
                    else:
                        new_node = mouseify_node(new_node)
                        converted[ ori_key ] = new_node
                    
                elif m2h : 
                    if ori_key in converted:
                        new_node = converted[ ori_key ]
                    else:
                        new_node = humanize_node(new_node)
                        converted[ ori_key ] = new_node
            
                i_keys[c] = new_node.key

                if ( new_node.key not in added_nKeys ):

                    if self._add_node( new_node ):
                        added_nKeys.add( new_node.key )
                    else:
                        node_was_filtered = True

            if not node_was_filtered:
                if data["interID"] not in added_iKeys :

                    thisinteraction = interaction( interID        =    data.get("interID"),
                                                   system         =    data.get("system"),
                                                   systemType     =    data.get("systemType"),
                                                   Author         =    data.get("Author"),
                                                   pmid           =    data.get("pmid"),
                                                   throughput     =    data.get("throughput"),
                                                   score          =    data.get("score"),
                                                   modification   =    data.get("modification"),
                                                   phenotypes     =    data.get("phenotypes"),
                                                   qualifications =    qualify if force_qualify else data.get('qualifications', qualify),
                                                   tags           =    data.get("tags"),
                                                   directed       =    directed,
                                                   srcDB          =    data.get("srcDB", ''),
                                                   nodeA          =    self.nodes[i_keys['A']],
                                                   nodeB          =    self.nodes[i_keys['B']])

                    if self._add_interaction( thisinteraction ):
                        added_iKeys.add(data["interID"])
                        
                else :
                    sys.stderr.write("\rWARNING: Interaction {} on line {} has already been added!\n".format(data["interID"],i)) ;

            dataline = infobj.readline() ;
            i       += 1;
            mu.waitbar( 80 * i / filelines, 80, showPct = True, fill = '%' )
            sys.stdout.flush()
            
        if (self.debug):
            sys.stderr.write("\nCompleted.\n")
        else : 
            mu.waitbar( 80, 80, fill = '%', showPct = True )
            sys.stderr.write('\n')


    def save(self, f, nodes = None, edges = None):

        if isinstance(f,str) :
            outfile = open(f,'w') ;
        else : 
            outfile = f ;


        if edges is None or type(edges) is not set : 
            edges = set(self.edges.values()) ;
        elif type(next(iter(edges))) is str : 
            edges = { self.edges[ek] for ek in edges if ek in self.edges }
        elif type(next(iter(edges))) == bgedge  : 
            pass ;

        if nodes is None or type(nodes) is not set : 
            nodes = set(self.nodes.values()) ;
        elif type(next(iter(nodes))) is str : 
            nodes = { self.nodes[nk] for nk in nodes if nk in self.nodes }
        elif isinstance(next(iter(nodes)),node) == node  : 
            pass ;

        for n in nodes :
            outfile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(n.key,n.entrez,
             n.biogrid,n.systematic,n.official,n.synonyms,n.organism)) ;

        outfile.write("\n") ; # this blank line will signal the difference between nodes and interaction

        interactors_outstring = "" ;
        for i in range(0,14) :
            interactors_outstring += "{}\t" ;
        interactors_outstring += "{}\n" ;

        #for i in list(self.the_data.values()) :
        for e in edges : 
            for i in e.interactions :
                if {i.nodeA,i.nodeB}.issubset(nodes) : 
                    dir_icon = '>' if i.directed else '.' ;
                    outfile.write(interactors_outstring.format(i.nodeA.key,i.nodeB.key,\
                     i.interID,i.system,i.systemType,i.Author,i.pmid,i.throughput,\
                     i.score,i.modification,i.phenotypes,i.qualifications,i.tags,\
                     i.srcDB,dir_icon)) ;
                else : 
                    pass ; 


    def load_from( self, infobj ) :
        """
            reading in a network file saved using the above 'save' method
            nodes and interactions are separated by an emtpy line
        """
        
        filelines   = sum( 1 for line in infobj ) ;
        infobj.seek(0) ; # return pointer to beginning of file

        i           = 0
        s           = infobj.readline()
        nodes       = True
        
        while s:
            if re.match( '^\s*$', s ):
                nodes = False
            else: 
                # reading nodes
                if nodes:
                    node_data = s.strip().split('\t') ;
                    if node_data[0] in self.nodes : 
                        if self.debug : 
                            sys.stderr.write('DEBUG> node'+node_data[0]+'_'+node_data[1]+' already exists in dataset.') ;
                        pass ; 
                    else : 
                        new_node  = node( key        = node_data[0],
                                          entrez     = node_data[1],
                                          biogrid    = node_data[2],
                                          systematic = node_data[3],
                                          official   = node_data[4],
                                          synonyms   = node_data[5],
                                          organism   = node_data[6])

                        self._add_node( new_node )

                # reading interactions
                else:
                    i_data        = s.strip().split('\t')

                    if self.nodes.get(i_data[0]) and self.nodes.get(i_data[1]): 

                        dir_state = True if i_data[14] == '>' else False
                        new_inter = interaction( nodeA          = self.nodes[i_data[0]],
                                                 nodeB          = self.nodes[i_data[1]],
                                                 interID        = i_data[2],
                                                 system         = i_data[3],
                                                 systemType     = i_data[4],
                                                 Author         = i_data[5],
                                                 pmid           = i_data[6],
                                                 throughput     = i_data[7],
                                                 score          = i_data[8],
                                                 modification   = i_data[9],
                                                 phenotypes     = i_data[10],
                                                 qualifications = i_data[11],
                                                 tags           = i_data[12],
                                                 srcDB          = i_data[13],
                                                 directed       = dir_state)
        
                        # add interaction to info for dataset
                        self._add_interaction( new_inter )
            
            i   += 1
            s    = infobj.readline()
            mu.waitbar( 80 * i / filelines, 80, showPct = True )

        infobj.close()

    def find_offs(self,official) : 

        out=list() ; 

        for n in self.nodes : 
            if n.official == official :
                out.append(n) ; 

        if len(out) == 0 : 
            return None ; 
        elif len(out) == 1 :
            return out.pop()
        else:
            return out ; 

    def clone(self,**kwargs) : 
        from io import StringIO 

        thebuffer=StringIO() ;
        self.save(thebuffer,nodes=kwargs.get('nodes')) ; 
        thebuffer.seek(0)
        other=dataSet(n_filter=kwargs.get('n_filter',self.n_filter),\
                      i_filter=kwargs.get('i_filter',self.i_filter)) ; 
        other.load_from(thebuffer) ; 

        return other ; 

    def __contains__(self,thing) : 

        if thing in self.nodes or thing in self.edges : 
            return True ;
        else : 
            return False ; 


        
def fetch_nodes(queries,thebiogrid) :

    node_queries=set() ; 

    for query in queries : # these are STRINGS, remember
        
        try : 
            node_queries.add(thebiogrid.nodes[query]) ;
        except KeyError : 
            pass ; 
            #sys.stderr.write("WARNING:    query {} has no node in database.\n".format(query)) ;

    return node_queries ;

def fetch_nodes_list(queries,thebiogrid) : 

    node_queries=list() ;

    for query in queries : 
        node_queries.append(thebiogrid.nodes.get(query,None)) ;

    return node_queries ; 
    

def fetch_nodes_regex(queries,thebiogrid,field="official") :

    qtype=type(next(iter(queries))) ; 
    regexset=set() ;

    if qtype is str :
        for query in queries : 
            regexset.add(re.compile(query,re.IGNORECASE)) ; 

    elif qytpe is REGEX_COMPILED_TYPE  :
        for query in queries : 
            regexset.add(query) ;
    else:
        sys.stderr.write("Container for queries does not hold items of type str or compiled regex.\n") ;
        raise ValueError

    node_queries=set() ;

    search_set=set(thebiogrid.nodes.values()) ;

    while search_set :
        s=search_set.pop() ;

        for r in regexset : 
            if r.match(s.__getattribute__(field)) :
                node_queries.add(s) ;
                break ;

    return node_queries ; 
            
def print_edges( edges, print_headers = True, print_mean_scores = False, print_weights = False, print_total_scores = False,
                 sep = '\t', fname = "", print_organisms = False, print_source = False, inter_string = 'pp', transform_scores = None,
                 print_quals = True, print_pps = False ):

    from numpy import log10

    if ( not fname ) :
        f=sys.stdout  ;
    #elif isfile(fname): 
        #sys.stdout.write('Appending to existing file {}\n'.format(fname)) ; 
        #f=open(fname,"a") ; 
    else : 
        f=open(fname,"w") ; 

    if ( print_headers) :
        f.write("Left{}Inter{}Right".format(sep,sep)) ;
        if ( print_quals ) :
            f.write("{}Qual".format(sep)) ;
        if ( print_weights ) :
            f.write("{}Wt".format(sep)) ;
        if ( print_mean_scores ) :
            f.write("{}Mean".format(sep)) ;
        if ( print_total_scores) :
            f.write("{}Total".format(sep)) ;
        if ( print_organisms ) :
            f.write("{}OrgA{}OrgB".format(sep,sep)) ; 
        if ( print_source ) :
            f.write("{}Source".format(sep,sep)) ; 
        if print_pps : 
            f.write("{}logp").format(sep) ; 

        f.write("\n") ;


    for edge in edges : 

        #edget=tuple(edge.nodes) ;
        if edge.directed : 
            thedir = '>' ; 
        else : 
            thedir = '^' ; 

        f.write("{}{}{}{}{}{}".format(edge.whence.official,sep,thedir,inter_string,sep,edge.to.official)) ;

        if ( print_quals ):
            f.write("{}{}".format(sep,edge.qual)) ;
        if ( print_weights ):
            f.write("{}{}".format(sep,edge.weight)) ;
        if ( print_mean_scores ):
            if transform_scores : 
                f.write("{}{}".format(sep,transform_scores(edge.meanscore))) ;
            else : 
                f.write("{}{}".format(sep,edge.meanscore)) ;
        if ( print_total_scores):
            if transform_scores : 
                f.write("{}{}".format(sep,transform_scores(edge.totalscore))) ;
            else : 
                f.write("{}{}".format(sep,edge.totalscore)) ;
        if ( print_organisms):
            if ( len(edge.nodes) == 2 ):
                f.write("{}{}{}{}".format(sep,edget[0].organism,sep,edget[1].organism)) ;
            else :
                f.write("{}{}{}{}".format(sep,edget[0].organism,sep,edget[0].organism)) ;
        if ( print_source):
            f.write("{}{}".format(sep,edge.source)) ;
        if print_pps : 
            f.write("{}{}".format(sep,-1*log10(edge.p)))

        f.write("\n") ;

    f.close()


def humanize_node( n ):    

    rbase.load( 'm2h' )
    rbase.load( 'hsg' )
    return convert_node( n, mouse_taxid, human_taxid, rbase.m2h, rbase.hsg )


def mouseify_node( n ):

    rbase.load( 'h2m' )
    rbase.load( 'mmg' )    
    return convert_node( n, human_taxid, mouse_taxid, rbase.h2m, rbase.mmg )
    

def convert_node(n, from_taxid, to_taxid, conv, targ ) :

    if n.organism == from_taxid : 

        if conv.get(n.entrez) is None and n.entrez in targ['EID']:
            target_entrez = n.entrez # the source entrezid is actually for the wrong org
        elif conv.get(n.entrez) is None:
            target_entrez = '' # source eid is not found
        elif len(conv[n.entrez]) > 1 :        
            for e in conv[n.entrez] :
                # try finding target with the same name
                if e in targ['EID']:
                    if targ['EID'][e]['Symbol'] == n.official.upper() : 
                        target_entrez = e
                        break ; 
            else :
                # pick target - the smallest entrezid                
                target_entrez = str(sorted(map(int, conv.get(n.entrez)))[0])

            sys.stderr.write('\rWARNING: Interactors.convert_node :'
                             +' {}:{} is not the unique homolog of {}:{}.\n'
                             .format( target_entrez, targ['EID'][target_entrez]['Symbol'],
                                      n.entrez, n.official ))

        else:
            # one to one mapping
            target_entrez = next(iter(conv[n.entrez])) ;

        if target_entrez :
            n.official = targ['EID'][target_entrez]['Symbol']
            n.organism = to_taxid
            n.entrez   = target_entrez
            n.key      = n.official+'_'+n.entrez

    return n ; 

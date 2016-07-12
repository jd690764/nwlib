import sys
import re
import time

from lib import config as cf

clue_list=["scarlet","mustard","green","white","peacock","plum"] ; 
clue_tuple=("scarlet","mustard","green","white","peacock","plum") ; 
clue_dict={ "miss"      : "scarlet",\
            "colonel"   : "mustard",\
            "mister"    : "green",\
            "maid"      : "white",\
            "missus"    : "peacock",\
            "professor" : "plum"} ;

clue_weapons=('knife','rope','lead_pipe','wrench','candlestick','revolver') ;
clue_rooms=('hall','lounge','ballroom','kitchen','dining_room','conservatory','libary','billiard_room','study') ;


def print_in_columns( thing, cols,colwidth=16,number=False,outfile=sys.stdout):
    i=0 ;
    for item in thing : 
        if ( number ):
            outs=" [{: >4}]{:->"+repr(colwidth)+"}" ; 
            outfile.write(outs.format(i,item)) ;
        else : 
            outs="{: >"+repr(colwidth)+"}"  ;
            outfile.write(outs.format(item)) ; 

        if ( i + 1 ) % cols  == 0 :
            outfile.write("\n") ; 
        i += 1;

    outfile.write("\n") ; 

def yesorno(prompt,yesses="Yy",noes="Nn") : 
    #yesnos=yesses+nos; 
    sys.stdout.write(prompt) ; 
    sys.stdout.flush() ;
    while True :
        userin=sys.stdin.readline() ; 

        if len(userin.strip()) == 0 :
            return True ; 
            break ; 
        else : 
            for c in userin : 
                if ( c in yesses ) :
                    return True ; 
                    break ; 
                elif ( c in noes ) :
                    return False ; 
                    break ; 
        sys.stderr.write("Invalid response.\n") ; 

def getLineCount(theFile):

        filelines=sum( 1 for line in theFile ) ;
        theFile.seek(0) ; # return pointer to beginning of file

        return filelines ; 

def waitbar(pos,end,showPct=False,fill='#') :
    pos=int(pos) ; 
    end=int(end) ; 
    outstr="\r|{:"+fill+"<"+str(pos+1)+"}{: >"+str(end-pos)+"}|" ;
    sys.stdout.write(outstr.format("","")) ;
    if ( showPct ) :
        sys.stdout.write(" ({:0>3}%)".format(100*pos//end)) ; 

    sys.stdout.flush() ;


def iterattr( collection,theAttr) :

    out=list() ;

    for i in collection :

        out.append(i.__getattribute__(theAttr)) ;

    return out ;

def stuff(theDict,theKey,theValue) : 

    if not theDict.get(theKey) : 
        theDict.update({ theKey : [ theValue ] }) ; 
        # add new item to field that does not yet exist
    elif type(theDict.get(theKey)) is list : 
        theDict[theKey].append(theValue);
        # append new gene to list already contained in self[field][newkey]
    else : 
        theDict.update({ theKey : [ theDict[theKey], theValue ] }) ;
        # init list with previous contents and with new gene


def stuffSet(theDict,theKey,theValue) : 

    if not theDict.get(theKey) : 
        theDict.update({ theKey : theValue }) ; 
        # add new item to field that does not yet exist
    elif type(theDict.get(theKey)) is set : 
        theDict[theKey].add(theValue);
        # append new gene to list already contained in self[field][newkey]
    else : 
        theDict.update({ theKey : { theDict[theKey], theValue } }) ;
        # init list with previous contents and with new gene

date6=lambda : time.strftime('%y%m%d') ;  

def randomProperNouns(numNouns) :

    from random import choice

    wordsf=file( cf.properNounsPath,'r') ;
    chosen_randoms=set() ;
    out=list() ;

    wordsl=wordsf.readlines() ; 

    while len(chosen_randoms) < numNouns :
        
        chosen_randoms.add(choice(wordsl)) ;

    for c in chosen_randoms : 

        out.append(c.strip()) ;

    wordsf.close() ;

    return out ; 

def indices(container) : 
    return range(0,len(container)) ;

def fields(*args,**kwargs) :

    sep=kwargs.get('sep','\t') ; 
    first=True ; 
    outstr='' ; 

    #print repr(sep) ; 

    for arg in args :

        if first : 
            first=False ; 
        else : 
            outstr += sep            

        if type(arg) is str : 
            outstr += arg ; 
        else:
            outstr += repr(arg) ; 
        

    if kwargs.get('withnewline',True) : 
        outstr += '\n' ; 

    return outstr ; 

def b4us(thestr) :

    outstr=''
    for c in thestr : 
        if c == '_' : 
            break ; 
        else : 
            outstr += c  ; 
    return outstr

def afterus(thestr) : 
    outstr='' ;
    getting=False
    for c in thestr : 
        if getting :
            outstr += c ; 
        elif c == '_' : 
            getting=True ;

    return outstr

def b4dot(thestr) : 

    outstr=''
    for c in thestr : 
        if c == '.' : 
            break ; 
        else : 
            outstr += c  ; 
    return outstr

def unstuff_values(stuffed_dict) :
    out=list() ; 
    for v in list(stuffed_dict.values()) : 
        if type(v) is list : 
            out.extend(v) ; 
        else :
            out.append(v) ; 

    return out

def unstuffset_values(stuffed_dict) : 
    out=set() ; 
    for v in list(stuffed_dict.values()) : 
        out.update(v) ; 

    return out

def absoverlap( symbols, taxid=9606 ): 

    from network.models import Entrez

    symbols = list(symbols) ; 
    entrez  = Entrez.objects.filter( symbol__in=symbols ).filter( taxid=taxid ).values( 'symbol', 'pubmed' )
    objs    = {}
    out     = dict()

    for entry in entrez:
        objs[ entry['symbol'] ] = set(entry['pubmed'].split(';'))

    for syma, pmidsa in objs.items() : 
        for symb, pmidsb in objs.items() : 
            if syma != symb:
                score = float(len( pmidsa & pmidsb )-1)**2 / float(len( pmidsa ^ pmidsb )  + 1 ) ; 
                out.update({ (syma,symb) : score }) ; 

    return out 

def peek(container) : 
    return next(iter(container)) ; 

def tabulate(fobj,sep='\t') :
    for line in fobj : 
        yield line.strip().split(sep) ; 

getty="Four score and seven years ago our fathers brought forth on this continent, a new nation, conceived in Liberty, and dedicated to the proposition that all men are created equal.\n\n Now we are engaged in a great civil war, testing whether that nation, or any nation so conceived and so dedicated, can long endure. We are met on a great battle-field of that war. We have come to dedicate a portion of that field, as a final resting place for those who here gave their lives that that nation might live. It is altogether fitting and proper that we should do this.\n\n But, in a larger sense, we can not dedicate -- we can not consecrate -- we can not hallow -- this ground. The brave men, living and dead, who struggled here, have consecrated it, far above our poor power to add or detract. The world will little note, nor long remember what we say here, but it can never forget what they did here. It is for us the living, rather, to be dedicated here to the unfinished work which they who fought here have thus far so nobly advanced. It is rather for us to be here dedicated to the great task remaining before us -- that from these honored dead we take increased devotion to that cause for which they gave the last full measure of devotion -- that we here highly resolve that these dead shall not have died in vain -- that this nation, under God, shall have a new birth of freedom -- and that government of the people, by the people, for the people, shall not perish from the earth."

def sortdict(thedict,reverse=False) : 

    sks=sorted(list(thedict.keys()),key=thedict.get,reverse=reverse) ; 

    out=[ (k,thedict[k]) for k in sks ]
   #out=list() ; 
   #for k in sks : 
   #    out.append((k,thedict.get(k)))

    return out

def getch():
    import termios
    import tty
    """getch() -> key character

    Read a single keypress from stdin and return the resulting character. 
    Nothing is echoed to the console. This call will block if a keypress 
    is not already available, but will not wait for Enter to be pressed. 

    If the pressed key was a modifier key, nothing will be detected; if
    it were a special function key, it may return the first character of
    of an escape sequence, leaving additional characters in the buffer.

    from http://code.activestate.com/recipes/577977-get-single-keypress/
    """
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(fd)
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)

    if ch == '\x03' : 
        raise KeyboardInterrupt

    return ch

import re


###
#find sequence
def find_seq_pattern(seqs, pattern=r'n.[stc]', ret_start=True):
    """
    find a pattern position in a sequence

    input
    ------
    seqs: string or list of strings to process

    pattern: regular expression to search

    ret_start: logical (default True)
        if true, return only the start of the pattern match position in pos
        if false, return (start,end) tuple of pattern match position in pos

    output
    ------
    pat: list of patterns matched in sequence
    pos: list of positions of patterns found

    if seqs is a list of sequences, then return is list of lists:
        pat[i][j]={pattern j found for sequence i}

    """
    #compile pattern
    p = re.compile(pattern)

    if type(seqs) is list:
        SL = seqs
    else:
        SL = [seqs]

    #containers for all matches
    patL = []
    posL = []

    #loop over all sequences
    for x in SL:
        pat = []
        pos = []
        #find all pattern matches
        for m in p.finditer(x):
            if m:
                pat.append(m.group())
                pos.append(m.span())

        if ret_start:
            pos = [x[0] for x in pos]

        patL.append(pat)
        posL.append(pos)

    if len(SL) == 1:
        patL = patL[0]
        posL = posL[0]
    return patL, posL


def pretty_seq(seq,
               pos,
               fmt=' %s ',
               colorcode=5,
               upper=True,
               lower=False,
               offset=3):
    """
    create a "pretty" sequence

    seq: string to consider
    
    pos: position to pretty up
        if pos[i] is tuple:
          pretty up pos[i][0]:pos[i][1]
        else:
          pretty up pos[i]:pos[i]+offset
    
    color: color code [0=black, 1=red, 2=green, 3=yellow, 4=blue, 5=magenta
                       6=cyan]
    
    offset: integer (default=3)
    
    """

    f = "\x1b[3%sm%s\x1b[0m" % (colorcode, fmt)

    a = ''
    x = seq
    start = 0

    for p in pos:
        if type(p) is tuple:
            ps, pe = p
        else:
            ps, pe = (p, p + offset)

            #not highlighted part
        a += x[start:ps]
        #highlighted part
        #a+="\x1b[35m %s \x1b[0m"%x[p:p+3]
        s = x[ps:pe]
        if upper:
            a += f % (s.upper())
        elif lower:
            a += f % (s.lower())
        else:
            a += f % s

        start = pe

    a += x[start:]
    return a

#http://ozzmaker.com/add-colour-to-text-in-python/

import colorama
import itertools

def seq_to_pretty_list(seq,
                       pos_list,
                       background_list,
                       fmt='%s',
                       text='BLACK',
                       upper=True,
                       lower=False,
                       offset=3):
    """
    create a "pretty" sequence

    seq: string to consider

    pos: position to pretty up
        if pos[i] is tuple:
          pretty up pos[i][0]:pos[i][1]
        else:
          pretty up pos[i]:pos[i]+offset

    text: string (default black)
        color in colorama (Note, converted to upper)

    backgroud : string (default YELLOW)

    offset: integer (default=3)

    """


    #f="\x1b[3%sm%s\x1b[0m"%(colorcode,fmt)
    #f="\033[1:40;43m%s\033[0m"%(fmt)
    #f = "\x1b[1;4%s;4%sm%s\x1b[0m" % (text, background, fmt)

    out = list(seq)
    for pos,color in itertools.zip_longest(pos_list,background_list):
        for p in pos:
            if type(p) is tuple:
                ps, pe = p
            else:
                ps, pe = (p, p+offset)
            #highlighted part
            for i in range(ps,pe):
                s=out[i]
                fmt = getattr(colorama.Fore, text.upper()) + getattr(colorama.Back, color.upper()) + '{}' + \
                    colorama.Style.RESET_ALL + colorama.Fore.BLACK
                x = fmt.format(s)
                
                out[i] = x
    
    return out


def extract_info_before_after(seq, pos, offset, offset_pos=3):
    """
    extract sequence info before and after pattern match
    
    input
    ------
    seq: sequence
    
    pos: position of matches (if list of ints, starts, if list if tuples, start/end)
    
    offset_pos: if pos is start, then offset of pattern match
    
    offset: offset for before/after strings
    
    output
    ------
    P: list of tuples of form (before,match,after)
    """
    P = []

    for p in pos:
        if type(p) is tuple:
            ps, pe = p
        else:
            ps, pe = (p, p + offset_pos)

        t = (seq[ps - offset:ps], seq[ps:pe], seq[pe:pe + offset])
        P.append(t)

    return P

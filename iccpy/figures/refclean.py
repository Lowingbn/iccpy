from sys import argv

def tex_bib_and_citations(tex_file):
    """ find all the citations in a tex file and find the bibliography """
    f = open(tex_file, 'r')
    bib_file = None
    bib_str = '\\bibliography{' 
    citations = []
    cite_str = 'cite'
    for line in f.readlines():
        if bib_str in line:
            idx1 = line.find(bib_str) + len(bib_str)
            idx2 = line[idx1:].find('}')
            bib_file = line[idx1:idx1+idx2]
        while cite_str in line:
            line = line[line.find(cite_str)+len(cite_str):]
            idx1 = line.find('{')
            idx2 = idx1+line[idx1:].find('}')
            this_cite = line[idx1+1:idx2].split(',')
            if len(this_cite)==0:
                raise Exception('bad citation!')
            citations.extend(c.strip() for c in this_cite)

    f.close()
    if bib_str is None:
        raise('No bibliography found!')
    bib_file = bib_file+'.bib'
    print 'bibliography:', bib_file
    return bib_file, citations

def bib_entries(bib_file):
    """ find all the entries in the bibliography """
    f = open(bib_file, 'r')
    bib = f.read()
    f.close()
    
    bib = bib.split('@')
    print 'Number of bibliography entries', len(bib)

    bib_all = {}
    for entr in bib:
        if '{' not in entr:
            continue
        idx1 = entr.find('{')
        idx2 = idx1+entr[idx1:].find(',')
        name = entr[idx1+1:idx2].strip()
        bib_all[name]='@'+entr
    return bib_all

def refclean(tex_file, out_file='cleanedbib.bib'):
    """ 
    make a new bibliography (out_file=cleaned.bib) 
    with only the references actually used in the paper
    """
    bib_file, citations = tex_bib_and_citations(tex_file)
    bib = bib_entries(bib_file)
    new_bib = {}
    print 'The following references are in', bib_file, 'but not cited in', tex_file,':'
    for name in bib.keys():
        if name not in citations:
            print name,
            continue
        new_bib[name]=bib[name]
    print '\n\nGenerating cleaned bibliography', out_file
    f = open(out_file, 'w')
    f.write('\n'.join(new_bib.values()))
    f.write('\n')
    f.close()
    print 'done.'

def bib_remove_duplicates(bib_file, out_file='cleanedbib.bib'):
    """ make a new bibliography (default cleanedbib.bib) with duplicates removed """
    all_bib = bib_entries(bib_file)
    f = open(out_file, 'w')
    f.write('\n'.join(all_bib.values()))
    f.write('\n')
    f.close()    

if __name__=='__main__':
    if len(argv)<2:
        print 'Usage refclean XYZ.tex|XYZ.bib'
        print 'For .tex files, searches for references and keeps only used ones in bibliography'
        print 'For .bib files, removes duplicates'
        exit(0)
        
    name = argv[1]
    if name[-4:]=='.tex':
        print 'latex file:', name
        refclean(name)
    elif name[-4:]=='.bib':
        print 'bibliography file:', name
        bib_remove_duplicates(name)
    else:
        print 'unkown file:', name
    

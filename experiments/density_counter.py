import sys
import truvari
import pysam

v = pysam.VariantFile(sys.argv[1])
m = truvari.Matcher()
m.params.pctseq = 0
m.params.sizemin = 20
m.params.sizefilt = 20
m.params.chunksize=200
c = truvari.chunker(m, ('base', v))

def get_bounds(cnk):
    mstart = sys.maxsize
    mend = 0
    for i in cnk:
        mstart = min(mstart, i.start)
        mend = max(mend, i.stop)
    return mstart, mend

for chunk, _ in c:
    s, e = get_bounds(chunk['base'])
    num = len(chunk['base'])
    print("%s\t%d\t%d\t%d" % (chunk['base'][0].chrom,  s, e, num))
    # Debugging for the sub-chunkers
    #from truvari.collapse import tree_size_chunker, tree_dist_chunker
    #for i, _ in tree_size_chunker(m, [(chunk, 0)]):
    #    if i['base']:
    #        s, e = get_bounds(i['base'])
    #        num = len(i['base'])
    #        print("%s\t%d\t%d\t%d\tsize" % (i['base'][0].chrom,  s, e, num))
    #    for j, _ in tree_dist_chunker(m, [(i, 0)]):
    #        if j['base']:
    #            s, e = get_bounds(j['base'])
    #            num = len(j['base'])
    #            print("%s\t%d\t%d\t%d\tdist" % (j['base'][0].chrom,  s, e, num))


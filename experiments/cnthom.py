import sys
import pysam

def cnthom(seq):
    maxspan = 5
    count = 0
    prev = None
    removed = 0
    for i in seq:
        if i == prev:
            count += 1
            if count >= maxspan:
                removed += 1
        else:
            count = 1
        prev = i
    return removed

v = pysam.VariantFile(sys.argv[1])
tot_bases = 0
hom_bases = 0
tot_var = 0
tot_any = 0
for entry in v:
    if entry.info['SVTYPE'] == 'DEL':
        seq = entry.ref
    else:
        seq = entry.alts[0]
    tot_bases += len(seq)
    ch = cnthom(seq)
    hom_bases += ch
    if ch:
        tot_any += 1
    tot_var += 1
print(hom_bases, '/', tot_bases, '=', hom_bases / tot_bases)
print(tot_any, '/', tot_var, '=', tot_any / tot_var)

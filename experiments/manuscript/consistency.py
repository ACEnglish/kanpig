"""
Given a multi-sample VCF with replicates, generate a report of consistent genotypes

Note, this currently has hard-coded FMR1 sample names.
"""
import pysam
import sys
import os
from collections import defaultdict, Counter
import truvari

bed = "/users/u233287/scratch/code/truvari/tickets/vd/GIAB-Q100.GRCh38.bed"
SIZEMIN=50
SIZEMAX=10000
VCFOUT = False

if not os.path.exists(sys.argv[1] + '.tbi'):
    sys.stderr.write("Index\n")
    exit(1)
vcf = pysam.VariantFile(sys.argv[1])
region_tree = truvari.build_region_tree(vcf, includebed=bed)
truvari.merge_region_tree_overlaps(region_tree)
vcf_i = truvari.region_filter(vcf, region_tree)

n_header = vcf.header.copy()
n_header.add_line('##INFO=<ID=NPAIR,Number=1,Type=Integer,Description="Number of present pair incosistent">')
if VCFOUT:
    out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)

samples = [_ for _ in list(vcf.header.samples) if _.startswith("FMR1")]
samples.sort()
lookup = defaultdict(list)
for i in samples:
    d = "-".join(i.split('-')[:2])
    lookup[d].append(i)

pairs = {k:v for k,v in lookup.items() if len(v) == 2}

counter = defaultdict(lambda: [0, 0, 0])
#nsamp = Counter()
for entry in vcf_i:
    if truvari.entry_variant_type(entry) not in [truvari.SV.DEL, truvari.SV.INS]:
        continue

    sz = truvari.entry_size(entry)
    if sz < SIZEMIN or sz > SIZEMAX:
        continue
    
    any_diff = False
    n_diff = 0
    for k, v in pairs.items():
        s1, s2 = v

        g1 = entry.samples[s1]['GT']
        s1_present = 1 in g1

        g2 = entry.samples[s2]['GT']
        s2_present = 1 in g2

        if not (s1_present or s2_present):
            continue

        identical = sum([_ for _ in g1 if _ is not None]) == sum([_ for _ in g2 if _ is not None])
        counter[k][0] += 1
        counter[k][1] += int(identical)
        counter[k][2] += int(s1_present and s2_present)
        if k == 'FMR1-2' and not identical:
            any_diff = True
        # This sample isn't the same individual
        if k != 'FMR1-1' and not identical:
            n_diff += 1
    #if any_diff:
    entry.translate(n_header)
    entry.info['NPAIR'] = n_diff
    if VCFOUT:
        out.write(entry)

    
if not VCFOUT:
    import pandas as pd
    d = pd.DataFrame(counter)
    d = d.T
    d.columns = ['num_variants', 'identical', 'both_present']
    print(d.to_csv(sep='\t'))
    #print(nsamp)

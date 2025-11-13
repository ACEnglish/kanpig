import sys
import truvari
import pandas as pd
import numpy as np

"""
Add parameters
"""
in_fn, out_fn = sys.argv[1:]
vcf = truvari.VariantFile(in_fn)
bed_fn = "/Users/english/code/SMaHT_MIMS/evaluate_tooling/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"

rows = []
for entry in vcf.fetch_bed(bed_fn):
    if entry.chrom in ['chrX', 'chrY'] \
        or entry.var_size() > 10000 \
        or entry.is_monrefstar() \
        or None in entry.samples[1]['GT']:
        continue

    b_gt = truvari.get_gt(entry.gt(0))
    o_gt = truvari.get_gt(entry.gt(1))

    rows.append([b_gt == o_gt, b_gt.name, o_gt.name, entry.samples[1]['DP'], *entry.samples[1]['AD'], entry.samples[1]['GQ']])
out = pd.DataFrame(rows, columns=['state', 'Ogt', 'Mgt', 'DP', 'AD_ref', 'AD_alt', 'GQ'])
out.to_csv(out_fn, index=False)

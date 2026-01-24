import truvari
import pandas as pd
import numpy as np
import sys

in_vcf, in_bed = sys.argv[1:]

vcf = truvari.VariantFile(in_vcf)

vaf_diff = []
counts = np.zeros((len(truvari.GT), len(truvari.GT)), dtype=int)
for entry in vcf.fetch_bed(in_bed):
    if entry.var_size() > 10000 \
            or None in entry.samples[0]['GT'] \
            or entry.is_monrefstar():
        continue
    vaf = entry.info['VAF'][1]
    germ_base_gt = truvari.get_gt(entry.samples['HG005']['GT'])
    germ_comp_gt = truvari.get_gt(entry.samples['SAMPLE']['GT'])
    counts[germ_base_gt.value, germ_comp_gt.value] += 1
    vaf_diff.append([germ_base_gt.name, vaf, entry.samples['SAMPLE']['AD'][1] / entry.samples['SAMPLE']['DP']])


same = np.diag(counts).sum()
check = np.sum(counts)
print(f"Overall Performance: {same} / {check} = {same / check * 100:.2f}%")

remove = counts[0, 0] # remove correct ref/ref to check germline performance
g_same = same - remove
g_check = check - remove
print(f"Germline Performance: {g_same} / {g_check} = {g_same / g_check * 100:.2f}%")

s_same = remove
s_check = np.sum(counts[0])
print(f"Somatic Performance: {s_same} / {s_check} = {s_same / s_check * 100:.2f}%")


labels = [_.name for _ in truvari.GT]
data = pd.DataFrame(counts, columns=labels, index=labels)
print(data)

vaf_diff = pd.DataFrame(vaf_diff, columns=['GT', 'base', 'comp'])
vaf_diff['Diff'] = vaf_diff['base'] - vaf_diff['comp']
print("Overall Expected - Predicted VAF")
print(vaf_diff['Diff'].describe())
print("By GT State Expected - Predicted VAF")
print(vaf_diff.groupby(['GT'])['Diff'].describe())


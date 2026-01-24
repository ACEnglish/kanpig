import truvari
import pandas as pd
import numpy as np
import sys
from collections import Counter

in_vcf, in_bed = sys.argv[1:]

vcf = truvari.VariantFile(in_vcf)

vaf_diff = []
samples = [f'samp{i}' for i in range(3)]
counts = np.zeros((len(truvari.GT), len(truvari.GT)), dtype=int)
per_counts = [np.zeros((len(truvari.GT), len(truvari.GT)), dtype=int) for _ in range(3)]

for entry in vcf.fetch_bed(in_bed):
    if entry.var_size() > 10000 \
            or None in entry.samples[0]['GT'] \
            or entry.is_monrefstar():
        continue
    vaf = entry.info['VAF'][1]
    germ_base_gt = truvari.get_gt(entry.samples['HG005']['GT'])
    tot_dp = 0
    tot_ad = 0
    seen_gts = Counter()
    for idx, i in enumerate(samples):
        if None in entry.samples[i]['GT']:
            continue
        m_gt = truvari.get_gt(entry.samples[i]['GT'])
        seen_gts[m_gt] += 1
        tot_dp += entry.samples[i]['DP']
        tot_ad += entry.samples[i]['AD'][1]
        per_counts[idx][germ_base_gt.value, m_gt.value] += 1

    
    if tot_dp:
        # I don't like just voting. Should be checking each one
        germ_comp_gt = seen_gts.most_common()[0][0]
        counts[germ_base_gt.value, germ_comp_gt.value] += 1
        vaf_diff.append([germ_base_gt.name, vaf, tot_ad / tot_dp])
    elif tot_ad == 0:
        germ_comp_gt = truvari.GT.UNK
        counts[germ_base_gt.value, germ_comp_gt.value] += 1
    else:
        germ_comp_gt = truvari.GT.NON
        counts[germ_base_gt.value, germ_comp_gt.value] += 1


same = np.diag(counts).sum()
check = np.sum(counts)
print("Cumulative Performance")
print(f"\tOverall: {same} / {check} = {same / check * 100:.2f}%")

remove = counts[0, 0] # remove correct ref/ref to check germline performance
g_same = same - remove
g_check = check - remove
print(f"\tGermline: {g_same} / {g_check} = {g_same / g_check * 100:.2f}%")

s_same = remove
s_check = np.sum(counts[0])
print(f"\tSomatic: {s_same} / {s_check} = {s_same / s_check * 100:.2f}%")

for idx, count in enumerate(per_counts):
    same = np.diag(count).sum()
    check = np.sum(count)
    print(f"samp{idx} Performance")
    print(f"\tOverall: {same} / {check} = {same / check * 100:.2f}%")

    remove = count[0, 0] # remove correct ref/ref to check germline performance
    g_same = same - remove
    g_check = check - remove
    print(f"\tGermline: {g_same} / {g_check} = {g_same / g_check * 100:.2f}%")

    s_same = remove
    s_check = np.sum(count[0])
    print(f"\tSomatic: {s_same} / {s_check} = {s_same / s_check * 100:.2f}%")


labels = [_.name for _ in truvari.GT]
data = pd.DataFrame(counts, columns=labels, index=labels)
print(data)

vaf_diff = pd.DataFrame(vaf_diff, columns=['GT', 'base', 'comp'])
vaf_diff['Diff'] = vaf_diff['base'] - vaf_diff['comp']
print("Cumulative Expected - Predicted VAF")
print(vaf_diff['Diff'].describe())
print("By GT State Expected - Predicted VAF")
print(vaf_diff.groupby(['GT'])['Diff'].describe())


"""
Given a VCF run through annotate.merge.py, generate a summary table of allele state performance metrics
"""
import pysam
import sys
import truvari
from collections import Counter

def print_report(tp, tn, fp, fn):
    compP = (tp + fp)
    baseP = (tp + fn)
    compN = (tn + fn)
    baseN = (tn + fp)
    ppv = tp / compP
    tpr = tp / baseP
    tnr = tn / baseN
    npv = tn / compN
    acc = (tp + tn) / (baseP + baseN)
    ba = (tpr + tnr) / 2
    f1 = 2 * ((ppv * tpr) / (ppv + tpr))
    print("compP", compP, sep='\t')
    print("baseP", baseP, sep='\t')
    print("compN", compN, sep='\t')
    print("baseN", baseN, sep='\t')
    print("ppv", ppv, sep='\t')
    print("tpr", tpr, sep='\t')
    print("tnr", tnr, sep='\t')
    print("npv", npv, sep='\t')
    print("acc", acc, sep='\t')
    print("ba", ba, sep='\t')
    print("f1", f1, sep='\t')

vcf = pysam.VariantFile(sys.argv[1])
cnt = Counter()

for entry in vcf:
    if truvari.entry_size(entry) > 10000:
        continue
    if 'ALState' not in entry.info:
        continue
    if entry.info["GTState"] not in ['Concordant', 'Discordant']:
        continue
    for i in entry.info["ALState"]:
        cnt[i] += 1

print_report(cnt["TP"], cnt["TN"], cnt["FP"], cnt["FN"])



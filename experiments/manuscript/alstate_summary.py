"""
Given a VCF run through annotate.merge.py, generate a summary table of allele state performance metrics
"""
import pysam
import sys


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
    print("compP", compP )
    print("baseP", baseP )
    print("compN", compN )
    print("baseN", baseN )
    print("ppv", ppv )
    print("tpr", tpr )
    print("tnr", tnr )
    print("npv", npv )
    print("acc", acc )
    print("ba", ba )
    print("f1", f1 )

vcf = pysam.VariantFile(sys.argv[1])

from collections import Counter

cnt = Counter()

for entry in vcf:
    if 'ALState' not in entry.info:
        continue
    for i in entry.info["ALState"]:
        cnt[i] += 1

print_report(cnt["TP"], cnt["TN"], cnt["FP"], cnt["FN"])



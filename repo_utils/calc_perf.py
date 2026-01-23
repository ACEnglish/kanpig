import truvari
import pandas as pd
import numpy as np
import sys

in_vcf, in_bed = sys.argv[1:]

vcf = truvari.VariantFile(in_vcf)

counts = np.zeros((len(truvari.GT), len(truvari.GT)), dtype=int)
for entry in vcf.fetch_bed(in_bed):
    if entry.var_size() > 10000 \
            or None in entry.samples[0]['GT'] \
            or entry.is_monrefstar(): # I don't like hard coding this
        continue
    a_gt = truvari.get_gt(entry.samples[0]['GT'])
    b_gt = truvari.get_gt(entry.samples[1]['GT'])
    counts[a_gt.value, b_gt.value] += 1


same = np.diag(counts).sum()
check = np.sum(counts)
print(f"Performance: {same} / {check} = {same / check * 100:.2f}%")
labels = [_.name for _ in truvari.GT]
data = pd.DataFrame(counts, columns=labels, index=labels)
print(data)

import sys
import pysam
from collections import Counter, defaultdict
import truvari
v = pysam.VariantFile(sys.argv[1])

stats = Counter()
gt_matrix = defaultdict(Counter)
for entry in v:
    if truvari.entry_is_filtered(entry) or truvari.entry_size(entry) < 20:
        continue
    stats['total'] += 1
    gt_s = entry.samples[0]['GT'] 
    gt_k = entry.samples[1]['GT'] 
    gt_matrix[gt_s][gt_k] += 1
    if entry.samples[0]['GT'] == entry.samples[1]['GT']:
        stats['exact'] += 1
    elif sum(entry.samples[0]['GT']) == sum(entry.samples[1]['GT']):
        stats['same_ac'] += 1
print(stats)
print('exact', stats['exact'] / stats['total'])
print('same_ac', stats['same_ac'] / stats['total'])
print('loose', (stats['exact'] + stats['same_ac']) / stats['total'])
import pandas as pd

x = {}
for key in gt_matrix:
    x[str(key)] = {}
    for key2 in gt_matrix[key]:
        x[str(key)][str(key2)] = gt_matrix[key][key2]
print(x)
print('gts')
print(pd.DataFrame(x))

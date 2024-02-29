"""
How good is cosine similarity at estimating sequence similarity?

I'm going to extract insertions from a VCF and if they're within 10bp, compare all vs all with sequence similarity and
with cosine similarity
"""
import kdp
import truvari
import pysam
import itertools
from scipy.spatial.distance import sqeuclidean
import numpy as np
def poor_jaccard(a, b):
    deno = 0
    neum = 0
    for i,j in zip(a, b):
        #if abs(i) + abs(j) < 2: continue
        deno += abs(i) + abs(j)
        neum += abs(i - j)
    if deno == 0:
        return 0
    if neum == 0:
        return 1
    return 1 - (neum / deno)

KMER=4
def work_group(cur_group):
    for i, j in itertools.combinations(cur_group, 2):
        sz1 = truvari.entry_size(i)
        sz2 = truvari.entry_size(j)
        seq_sim = truvari.seqsim(i.alts[0], j.alts[0])
        seq_sim2 = truvari.entry_seq_similarity(i, j)
        szsim = truvari.entry_size_similarity(i, j)
        k1 = kdp.Haplotype.from_vcf(i, KMER)
        k2 = kdp.Haplotype.from_vcf(j, KMER)
        cos_sim = kdp.cosinesim(k1.kfeat, k2.kfeat)
        wcos_sim = kdp.weighted_cosinesim(k1.kfeat, k2.kfeat)
        esim = sqeuclidean(k1.kfeat, k2.kfeat)
        pj = poor_jaccard(k1.kfeat, k2.kfeat)
        print(i.chrom, i.pos, j.pos, round(sz1, 4), 
                                     round(sz2, 4),
                                    round(seq_sim,4),
                                    round(seq_sim2,4),
                                    round(cos_sim,4),
                                    round(wcos_sim,4),
                                    round(szsim[0],4),
                                    round(esim,4),
                                    round(pj, 4),
                                    sep='\t')


fn = "/Users/english/code/kfdphase/kdprs/test/test2.vcf.gz"
fn = "/Users/english/code/truvari/tickets/vcfdist_reeval/Q100-dipcall.HG002.GRCh38.vcf.gz"
v = pysam.VariantFile(fn)
print("chrom\tpos1\tpos2\tsz1\tsz2\tseqsim\tunroll\tcossim\twcossim\tszsim\teucsim\tpj")
last_pos = 0
cur_group = []
for entry in v:
    if truvari.entry_variant_type(entry) != truvari.SV.INS:
        continue
    if truvari.entry_size(entry) < 20:
        continue
    if entry.pos - last_pos > 10:
        work_group(cur_group)
        cur_group = []
    last_pos = entry.pos
    cur_group.append(entry)

work_group(cur_group)

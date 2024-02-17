"""
Methods for creating haplotypes from either a phased base vcf or long reads

TODO:
    Need to standardize the return type
    Then update the kdp_job so that it'll make the struct
    Then update phase_region so that it'll consume the struct
    Then I just have to make the bam parser do the same thing
"""
from collections import defaultdict
import pysam

def vcf_haps(variants, kmer):
    """
    Parse a set of phased variants
    """
    h1_k = kdp.seq_to_kmer("", kmer)
    h2_k = kdp.seq_to_kmer("", kmer)
    h1_s = 0
    h2_s = 0
    h1_n = 0
    h2_n = 0
    for v in variants:
        k, sz = kdp.var_to_kfeat(v, kmer)
        if v.samples[0]['GT'][0] == 1:
            h1_k += k
            h1_s += sz
            h1_n += 1
        if len(v.samples[0]['GT']) > 1 and v.samples[0]['GT'][1] == 1:
            h2_k += k
            h2_s += sz
            h2_n += 1
    return h1_k, h1_s, h1_n, h2_k, h2_s, h2_n


def bam_haps(bam, chrom, start, end):
    """
    Pileup a region and return same thing as vcf_haps
    """
    #chrom, start, end = "chr20", 20827970, 20827980 # insertion
    all_cov = 0
    for pileup_column in bam.pileup(chrom, start, end, truncate=True):
        # Check for deletions
        all_cov += pileup_column.n
        for pileup_read in pileup_column.pileups:
            if pileup_read.indel > 20:  # Insertion greater than 20 bp
                print(pileup_column.reference_pos, pileup_read.indel)
                print(pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileup_read.indel])
            elif pileup_read.indel < -20:  # Deletion greater than 20 bp
                print(pileup_column.reference_pos, pileup_read.indel)

    print(all_cov / (end - start))

    

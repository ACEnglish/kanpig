"""
Methods for creating haplotypes from either a phased base vcf or long reads

TODO:
    Need to standardize the return type
    Then update the kdp_job so that it'll make the struct
    Then update phase_region so that it'll consume the struct
    Then I just have to make the bam parser do the same thing
"""
import copy
from dataclasses import dataclass
from collections import defaultdict

import pysam
from numpy.typing import ArrayLike

import kdp

@dataclass
class Haplotype():
    """
    Holds the kfeat and size deltas of a sequence
    """
    kfeat: ArrayLike # kmer featurization of this haplotype
    size: int # Size delta of this haplotype
    n: int # Number of changes in this haplotype
    
    @staticmethod
    def new(kmer):
        return Haplotype(kdp.seq_to_kmer("", kmer), 0, 0)

    @staticmethod
    def from_vcf(entry, kmer=3):
        """
        Turn variant record into a kfeat
        """
        alt = kdp.seq_to_kmer(entry.alts[0], kmer)
        ref = kdp.seq_to_kmer(entry.ref, kmer)
        szdiff = len(entry.alts[0]) - len(entry.ref)
        return Haplotype(alt - ref, szdiff, 1)

    def __iadd__(self, other):
        """
        self += other modifies self
        """
        self.kfeat += other.kfeat
        self.size += other.size
        self.n += other.n
        return self

    def __add__(self, other):
        """
        new = self + other returns copy of self
        """
        ret = copy.deepcopy(self)
        ret += other
        return ret

def vcf_haps(variants, kmer):
    """
    Parse a set of phased variants and return the two Haplotypes
    """
    h1 = Haplotype.new(kmer)
    h2 = Haplotype.new(kmer)
    for entry in variants:
        m_hap = Haplotype.from_vcf(entry, kmer)
        if entry.samples[0]['GT'][0] == 1:
            h1 += m_hap
        if len(entry.samples[0]['GT']) > 1 and entry.samples[0]['GT'][1] == 1:
            h2 += m_hap
    return h1, h2


def bam_haps(bam, chrom, start, end):
    """
    Pileup a region and return same thing as vcf_haps
    """
    # chrom, start, end = "chr20", 20827970, 20827980 # insertion
    all_cov = 0
    for column in bam.pileup(chrom, start, end, truncate=True):
        # Check for deletions
        all_cov += column.n
        for read in column.pileups:
            if read.indel > 20:  # Insertion greater than 20 bp
                print(column.reference_pos, read.indel)
                print(read.alignment.query_sequence[read.query_position:read.query_position + read.indel])
            elif read.indel < -20:  # Deletion greater than 20 bp
                print(column.reference_pos, read.indel)
                # Need to fetch seq from reference, probably going to send reference in with bam to save IO
    # Region coverage
    print(all_cov / (end - start))

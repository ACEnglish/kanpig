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
from sklearn.cluster import KMeans
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
    coverage: int = 1 # How many sequences support this haplotype

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
        Adding is only for extending a single sequence
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


def bam_haps(bam, chrom, start, end, kmer=4, buffer=100, sizemin=20):
    """
    Pileup a region and return same thing as vcf_haps
    """
    tot_cov = 0
    m_haps = {} # readname: Haplotype
    for column in bam.pileup(chrom, start - buffer, end + buffer, truncate=True):
        # Check for deletions
        for read in column.pileups:
            # Guard against partial alignments which mess up the kfeat 
            # Will revisit when I can turn a Haplotype into a single-path graph
            if not ((read.alignment.reference_start < start) and (read.alignment.reference_end > end)):
                continue
            # Only count reads that span the region towards coverage
            tot_cov += 1
            # Only consider things greater than 20bp
            if abs(read.indel) < sizemin:
                continue

            if (read.indel ^ 1) > 0:  # Insertion
                seq = read.alignment.query_sequence[read.query_position:read.query_position + read.indel]
                hap = kdp.Haplotype(kdp.seq_to_kmer(seq, kmer), read.indel, 1)
            else:  # Deletion
                m_start = column.reference_pos - start
                m_end = m_start + abs(read.indel)
                hap = kdp.Haplotype(-kdp.seq_to_kmer(refseq[m_start: m_end], kmer), read.indel, 1)

            if read.alignment.query_name not in m_haps:
                m_haps[read.alignment.query_name] = hap
            else:
                m_haps[read.alignment.query_name] += hap

    # Region coverage
    coverage = tot_cov / (end - start)
    hap1, hap2 = read_cluster(m_haps, kmer, coverage)
    return hap1, hap2

def average_of_haps(m_haps, kmer):
    """
    Make a new hap that's the average of all the reads
    """
    ks, sz, n = list(zip(*[_.kfeat, _.size, _.n for _ in m_haps]))
    hap = kdp.Haplotype.new(kmer)
    hap.kfeat = np.mean(ks, axis=0)
    hap.size = np.mean(sz)
    hap.n = np.mean(n)
    hap.coverage = len(m_haps)
    return hap

def read_cluster(m_haps, kmer, coverage):
    """
    if the number of variant reads is appx 0%:
        we expect reference homozygous

    if number of variant reads is approximately 50%:
        we expect only 1 alternate group
        and we'll return a reference allele (as hap1)

    if the number of variant reads is like ≥80%:
        we expect either a compound het, or a hom
        We need to 
    if the number of variant reads is approximately 0%, we expect reference homozygous across
    # If we 
    Todo - The coverage from non-variant reads is the reference haplotype
    """
    # Add in the reference for clustering... 
    # Assuming no read is named reference...
    alt_proportion = len(all_ks) / coverage
    if alt_proportion < 0.85: # HARD CODED - we need to allow one of the paths to be reference
        all_ks['reference'] = kdp.Haplotype(kdp.seq_to_kmer("", 4), 0, reg_cov - len(all_ks))

    # clustering
    kmeans = KMeans(n_clusters=2, random_state=0)
    weight = [_.coverage for _ in all_ks.values()]
    # Separate reads into haplotypes
    grps = list(zip(kmeans.fit_predict([_.kfeat for _ in all_ks.values()], sample_weight=weight), all_ks.keys()))
    ref_zero = None
    haps = [[], []]
    for m_gt, hap_name in grps:
        if hap_name == 'reference':
            ref_is_zero = m_gt == 0
        haps[m_gt].append(all_ks[hap_name])

    if ref_is_zero is None: # possibly a compound het
        hap1 = average_of_haps(haps[0], kmer)
        hap2 = average_of_haps(haps[0], kmer)
    elif ref_is_zero:
        hap1 = all_ks['reference']
        # things that went into haps[0] are now lost...
        hap2 = average_of_haps(haps[1], kmer)
    else:
        hap1 = all_ks['reference']
        # things that went into haps[1] are now lost...
        hap2 = average_of_haps(haps[0], kmer)
    
    return hap1, hap2

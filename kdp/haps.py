"""
Methods for creating haplotypes from either a phased base vcf or long reads

TODO:
    Need to standardize the return type
    Then update the kdp_job so that it'll make the struct
    Then update phase_region so that it'll consume the struct
    Then I just have to make the bam parser do the same thing
"""
import math
import copy
import logging
from dataclasses import dataclass
from collections import defaultdict

import pysam
import numpy as np
from numpy.typing import ArrayLike
from sklearn.cluster import KMeans
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning

import kdp

# If non-ref sequence takes less than this percent of the coverage
# we are expecting a reference haplotype
REFTHRESHOLD=0.85

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
    def new(kmer, coverage=1):
        return Haplotype(kdp.seq_to_kmer("", kmer), 0, 0, coverage)

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

    def __eq__(self, other):
        return np.all(self.kfeat == other.kfeat)

    def __hash__(self):
        return hash(str(self.kfeat))

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


def bam_haps(bam, refseq, chrom, reg_start, reg_end, kmer=4, buffer=100, cossim=0.95, sizemin=20, sizemax=50000):
    """
    Pileup a region and return same thing as vcf_haps
    """
    tot_cov = 0
    # TODO - I don't like this buffer business.. maybe
    m_haps = {} # readname: Haplotype
    for column in bam.pileup(chrom, reg_start - buffer, reg_end + buffer, truncate=True):
        # Check for deletions
        tot_cov += column.n
        for read in column.pileups:
            # Guard against partial alignments which mess up the kfeat 
            # Will revisit when I can turn a Haplotype into a single-path graph
            if not ((read.alignment.reference_start < reg_start) and (read.alignment.reference_end > reg_end)):
                continue
            # Only count reads that span the region towards coverage
            # Only consider things greater than 20bp
            if not (sizemin <= abs(read.indel) <= sizemax):
                continue

            if (read.indel ^ 1) > 0:  # Insertion
                seq = read.alignment.query_sequence[read.query_position:read.query_position + read.indel]
                hap = kdp.Haplotype(kdp.seq_to_kmer(seq, kmer), read.indel, 1)
                logging.debug('INS %d @ %d -> %s', len(seq), read.query_position, seq)
            else:  # Deletion
                # Add 1 for the .. reason..
                m_start = column.reference_pos - (reg_start - buffer) + 1
                m_end = m_start + abs(read.indel)
                seq = refseq[m_start: m_end]
                hap = kdp.Haplotype(-kdp.seq_to_kmer(seq, kmer), read.indel, 1)
                logging.debug('DEL %d @ %d -> %s', len(seq), column.reference_pos, seq)

            if read.alignment.query_name not in m_haps:
                m_haps[read.alignment.query_name] = hap
            else:
                m_haps[read.alignment.query_name] += hap
    # Region coverage
    coverage = int(tot_cov / (reg_end - reg_start + 2*buffer))
    logging.debug('reg coverage %d', coverage)
    # Nothing to compare
    if coverage == 0 or len(m_haps) == 0:
        ret = kdp.Haplotype.new(kmer, coverage)
        return ret, ret

    # No alternates
    #if len(m_haps) == 0:
        #return [kdp.Haplotype.new(kmer, coverage)]

    m_haps = hap_deduplicate(m_haps.values())

    # Old clustering... which had its benefits, but... like... yeah
    # no reads or no alts means reference
    # TODO: no coverage means we can't genotype it. Separate out these two conditions
    if (coverage == 0) or (len(m_haps) == 0):
        all_ref = kdp.Haplotype.new(kmer, coverage)
        return all_ref, all_ref
    elif len(m_haps) == 1: # one read can't be clustered, just return it
        hap2 = list(m_haps.values())[0]
        if (hap2.coverage / coverage) < REFTHRESHOLD:
            hap1 = kdp.Haplotype.new(kmer, coverage)
        else:
            hap1 = hap2
        return hap1, hap2

    hap1, hap2 = read_cluster(m_haps, kmer, coverage, cossim)
    return hap1, hap2

def hap_deduplicate(m_haps):
    """
    """
    ret = {}
    for i in m_haps:
        if i in ret:
            ret[i].coverage += 1
        else:
            ret[i] = i
    return ret

def average_of_haps(m_haps):
    """
    Return the hap from a cluster of haps the highest coverage or fewest changes
    #Make a new hap that's the average of all the reads
    """
    # Pick the one with the highest coverage or closest to the centroid
    best = m_haps[0]
    for i in m_haps[1:]:
        if i.coverage > best.coverage:
            best = i
        elif i.coverage == best.coverage and i.n < best.n:
            best = i
    for i in m_haps:
        if i != best:
            best.coverage += i.coverage
    return best

@ignore_warnings(category=ConvergenceWarning)
def read_cluster(all_ks, kmer, coverage, cossim):
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
    # clustering
    kmeans = KMeans(n_clusters=2, random_state=0)
    weight = [_.coverage for _ in all_ks.values()]
    # Separate reads into haplotypes
    grps = kmeans.fit_predict([_.kfeat for _ in all_ks.values()], sample_weight=weight)
    n_grps = len(set(grps))
    alt_cov = sum([_.coverage for _ in all_ks])
    # REF & HET
    if n_grps == 1 and (alt_cov / coverage) < REFTHRESHOLD :
        logging.debug('ref_het')
        hap1 = kdp.Haplotype(kdp.seq_to_kmer("", 4), 0, coverage - alt_cov)
        hap2 = average_of_haps(all_ks.values())
        return hap1, hap2

    # HOM ALT
    if n_grps == 1:
        logging.debug('hom_alt')
        hap = average_of_haps(all_ks.values())
        return hap, hap

    # Compound HET?
    logging.debug('compound_het')
    haps = [[], []]
    for m_gt, hap_name in zip(grps, all_ks.keys()):
        haps[m_gt].append(all_ks[hap_name])

    hap1 = average_of_haps(haps[0])
    hap2 = average_of_haps(haps[1])

    # if the two haps are super similar, just combine them
    # and then double check the reference threshold
    if kdp.cosinesim(hap1.kfeat, hap2.kfeat) >= cossim:
        hap2 = average_of_haps([hap1, hap2])
        if (alt_cov / coverage) < REFTHRESHOLD:
            hap1 = kdp.Haplotype(kdp.seq_to_kmer("", 4), 0, coverage - alt_cov)
        else:
            hap1 = hap2
        return hap1, hap2

    # Keep Compound Het
    return hap1, hap2


def genotyper(totCov, altCov, priors=None):
    """
    Given total coverage and altCoverage, try to calculate how many copies
    of the alt are at the position (ref/het/hom) and a quality score
    returns two lists
    - probabilities of 0, 1, or two copies of the allele at the location
    - phred-scaled quality scores of those probs
    """
    #We have no information.. should give up
    if totCov == 0:
        return None

    # previously had avgCov
    if priors is None:
        priors = [0.05, 0.5, 0.95]

    # if len(priors) != 3: # raise exception?

    def log_choose(n, k):
        """ swap for efficiency if k is more than half of n """
        r = 0.0
        if k * 2 > n:
            k = n - k

        for d in range(1, k + 1):
            r += math.log(n, 10)
            r -= math.log(d, 10)
            n -= 1

        return r

    total = totCov  # refCoverage + altCoverage if avgCov is None else avgCov
    alt = altCov  # int(spot.tags["szCount"])
    non_alt = total - alt

    gtList = []

    comb = log_choose(total, alt)
    for p_alt in priors:
        gtList.append(comb + alt * math.log(p_alt, 10) + non_alt * math.log(1 - p_alt, 10))

    return gtList

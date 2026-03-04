use crate::kplib::{pileup::ReadPileup, vcftraits::Svtype, KmerVec, SequenceMeta};
use itertools::Itertools;
use std::{
    cmp::Ordering,
    fmt::{Debug, Formatter, Result},
    hash::{Hash, Hasher},
};

#[derive(Clone)]
pub struct Haplotype {
    pub size: i64,
    pub n: u64,
    pub kmers: KmerVec,
    pub parts: Vec<(i64, KmerVec)>,
    pub partial: usize,
    pub meta: SequenceMeta,
}

impl Haplotype {
    pub fn new(kmers: KmerVec, size: i64, n: u64, hap_meta: SequenceMeta) -> Self {
        Self {
            size,
            n,
            kmers: kmers.clone(),
            parts: vec![(size, kmers)],
            partial: 0,
            meta: hap_meta,
        }
    }

    // Create an empty haplotype
    pub fn blank(meta: SequenceMeta) -> Haplotype {
        Haplotype {
            size: 0,
            n: 0,
            kmers: KmerVec::blank(),
            parts: vec![],
            partial: 0,
            meta,
        }
    }

    pub fn from_readpileup(pileup: ReadPileup, kmer: (u8, u8)) -> Haplotype {
        let mut ret = Haplotype::blank(pileup.meta.clone());
        for p in pileup.pileups.iter() {
            ret.size += p.size;
            ret.n += 1;

            let other_kfeat = KmerVec::new(
                p.sequence
                    .as_ref()
                    .expect("You didn't fill in pileup sequence"),
                kmer,
                p.indel == Svtype::Del,
            );

            ret.kmers += &other_kfeat;

            ret.parts.push((p.size, other_kfeat));
        }
        ret
    }

    /// Clear the Metadata and return a clone
    /// This allows us to preserve the KmerVec but update the meta
    pub fn clear_clone(&self, id: usize) -> Haplotype {
        let mut ret = self.clone();
        let mut n_meta = SequenceMeta::new_blank(self.meta.coverage.len());
        n_meta.id = id;
        ret.meta = n_meta;
        ret
    }

    /// Add another variant to a Haplotype
    /// This is useful for combining "haplotypes" that are actually sub-haplotypes
    /// e.g. variants across a read
    pub fn add(&mut self, other: &Haplotype) {
        self.kmers += &other.kmers;
        self.size += other.size;
        self.n += 1;
        self.parts.push((other.size, other.kmers.clone()));
    }

    /// Create new haplotypes of subsets of the variants
    /// This is essentially allowing for false negatives in the graph by pretending
    /// the haplotype doesn't have all the variants, and if so, perhaps there is a
    /// better fit.
    pub fn partial_haplotypes(&self, max_fns: usize, max_parts: usize) -> Vec<Haplotype> {
        let mut ret = vec![];
        let m_len = self.parts.len();
        if m_len >= max_parts {
            ret.push(self.clone());
            return ret;
        }
        let lower = if m_len <= max_fns { 1 } else { m_len - max_fns };
        for n in (lower..(m_len + 1)).rev() {
            for subset in self.parts.iter().combinations(n) {
                let mut cur_hap = Haplotype::blank(self.meta.clone());
                for part in subset.iter() {
                    cur_hap.size += part.0;
                    cur_hap.kmers += &part.1;
                    cur_hap.n += 1;
                }
                cur_hap.partial = m_len - n;
                ret.push(cur_hap);
            }
        }
        ret
    }
}

impl PartialOrd for Haplotype {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Haplotype {
    fn cmp(&self, other: &Self) -> Ordering {
        let coverage_ordering = self
            .meta
            .coverage
            .iter()
            .sum::<u64>()
            .cmp(&other.meta.coverage.iter().sum::<u64>());
        if coverage_ordering != Ordering::Equal {
            return coverage_ordering;
        }

        let changes_ordering = self.n.cmp(&other.n).reverse();
        if changes_ordering != Ordering::Equal {
            return changes_ordering;
        }

        let size_ordering = self.size.cmp(&other.size);
        if size_ordering != Ordering::Equal {
            return size_ordering;
        }

        self.kmers.cmp(&other.kmers)
    }
}

impl PartialEq for Haplotype {
    fn eq(&self, other: &Self) -> bool {
        self.meta.coverage.iter().sum::<u64>() == other.meta.coverage.iter().sum::<u64>()
            && self.size == other.size
            && self.n == other.n
            && self.kmers == other.kmers
    }
}

impl Eq for Haplotype {}

impl Hash for Haplotype {
    // Maybe should take other Haplotype properties into account?
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.kmers.hash(state);
    }
}

impl Debug for Haplotype {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.debug_struct("Haplotype")
            .field("size", &self.size)
            .field("n", &self.n)
            .field("coverage", &self.meta.coverage)
            .field("ps", &self.meta.ps)
            .field("hp", &self.meta.hp)
            .field("samp", &self.meta.samples_flag)
            // Exclude kmers from the debug output
            .finish()
    }
}

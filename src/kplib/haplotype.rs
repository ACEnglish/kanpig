use crate::kplib::KmerVec;
use itertools::Itertools;
use std::{
    cmp::Ordering,
    fmt::{Debug, Formatter, Result},
    hash::{Hash, Hasher},
};

/// This holds read information that's eventually passed to PathScore
/// And then is used by the GenotypeAnno to fill in FORMAT fields
/// In order to allow reads across samples to talk to one another
/// we need to use Vectors. For a single sample operation, we will be
/// accessing everything simply as attribute[0]. For multi-sample, we'll
/// use e.g. attribute[0] for proband, attribute[1] for mother, etc.
/// Since reads can consolidate into a single haplotype, we use the
/// sample_flag as a shortcut to know what samples contributed to the
/// haplotype. e.g. flag & 1 means this is a proband haplotype
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct HaplotypeMeta {
    // Once we cluster/genotype, we need to be able to identify them
    // esp for mosaic
    pub id: usize,
    pub coverage: Vec<u64>,
    pub ps: Vec<Option<u32>>,
    pub hp: Vec<Option<u8>>,
    pub samples_flag: usize,
    pub rnames: Vec<String>,
}

impl HaplotypeMeta {
    pub fn new(sample_idx: usize, num_samples: usize) -> Self {
        let mut coverage = vec![0u64; num_samples];
        coverage[sample_idx] += 1;
        let rnames = vec![];

        HaplotypeMeta {
            id: 0,
            coverage,
            ps: vec![None; num_samples],
            hp: vec![None; num_samples],
            samples_flag: 2_usize.pow(sample_idx as u32),
            rnames,
        }
    }

    /// New HaplotypeMetadata that doesn't belong to anyone
    pub fn new_blank(num_samples: usize) -> Self {
        let coverage = vec![0u64; num_samples];
        HaplotypeMeta {
            id: 0,
            coverage,
            ps: vec![None; num_samples],
            hp: vec![None; num_samples],
            samples_flag: 0,
            rnames: vec![],
        }
    }

    pub fn combine(&mut self, other: &HaplotypeMeta) {
        for (self_cov, other_cov) in self.coverage.iter_mut().zip(&other.coverage) {
            *self_cov += other_cov;
        }

        for (self_ps, other_ps) in self.ps.iter_mut().zip(&other.ps) {
            if self_ps.is_none() {
                *self_ps = *other_ps;
            }
        }

        for (self_hp, other_hp) in self.hp.iter_mut().zip(&other.hp) {
            if self_hp.is_none() {
                *self_hp = *other_hp;
            }
        }

        self.rnames.extend(other.rnames.clone());

        self.samples_flag |= other.samples_flag;
    }
}

#[derive(Clone)]
pub struct Haplotype {
    pub size: i64,
    pub n: u64,
    pub kmers: KmerVec,
    pub parts: Vec<(i64, KmerVec)>,
    pub partial: usize,
    pub meta: HaplotypeMeta,
}

impl Haplotype {
    pub fn new(kmers: KmerVec, size: i64, n: u64, hap_meta: HaplotypeMeta) -> Self {
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
    pub fn blank(meta: HaplotypeMeta) -> Haplotype {
        Haplotype {
            size: 0,
            n: 0,
            kmers: KmerVec::blank(),
            parts: vec![],
            partial: 0,
            meta,
        }
    }

    /// Clear the Metadata and return a clone
    /// This allows us to preserve the KmerVec but update the meta
    pub fn clear_clone(&self, id: usize) -> Haplotype {
        let mut ret = self.clone();
        let mut n_meta = HaplotypeMeta::new_blank(self.meta.coverage.len());
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

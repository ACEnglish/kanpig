use crate::kplib::seq_to_kmer;
use itertools::Itertools;
use std::{
    cmp::Ordering,
    fmt::{Debug, Formatter, Result},
    hash::{Hash, Hasher},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct HaplotypeMeta {
    pub coverage: u64,
    pub ps: Option<u32>,
    pub hp: Option<u8>,
    pub samples_idx: u8,
}
impl HaplotypeMeta {
    pub fn new(samples_idx: u8) -> Self {
        HaplotypeMeta {
            coverage: 1,
            ps: None,
            hp: None,
            samples_idx,
        }
    }
}

#[derive(Clone)]
pub struct Haplotype {
    pub size: i64,
    pub n: u64,
    pub kfeat: Vec<f32>,
    pub parts: Vec<(i64, Vec<f32>)>,
    pub partial: usize,
    pub meta: HaplotypeMeta,
}

impl Haplotype {
    pub fn new(kfeat: Vec<f32>, size: i64, n: u64, hap_meta: HaplotypeMeta) -> Self {
        Self {
            size,
            n,
            kfeat: kfeat.clone(),
            parts: vec![(size, kfeat)],
            partial: 0,
            meta: hap_meta,
        }
    }

    // Create an empty haplotype
    pub fn blank(kmer: u8, meta: HaplotypeMeta) -> Haplotype {
        let mk = seq_to_kmer(&[], kmer, false, 0);
        Haplotype {
            size: 0,
            n: 0,
            kfeat: mk.clone(),
            parts: vec![],
            partial: 0,
            meta,
        }
    }

    // Add another variant to a Haplotype
    pub fn add(&mut self, other: &Haplotype) {
        if !self.kfeat.len() == other.kfeat.len() {
            panic!("Cannot add haplotypes of different kmer size");
        }
        self.kfeat
            .iter_mut()
            .zip(other.kfeat.iter())
            .for_each(|(x, y)| *x += y);
        self.size += other.size;
        self.n += 1;
        self.parts.push((other.size, other.kfeat.clone()));
    }

    pub fn partial_haplotypes(&self, kmer: u8, max_fns: usize, max_parts: usize) -> Vec<Haplotype> {
        let mut ret = vec![];
        let m_len = self.parts.len();
        if m_len >= max_parts {
            ret.push(self.clone());
            return ret;
        }
        let lower = if m_len <= max_fns { 1 } else { m_len - max_fns };
        for i in (lower..(m_len + 1)).rev() {
            for j in self.parts.iter().combinations(i) {
                let mut cur_hap = Haplotype::blank(kmer, self.meta);
                for k in j.iter() {
                    cur_hap.size += k.0;
                    cur_hap
                        .kfeat
                        .iter_mut()
                        .zip(k.1.iter())
                        .for_each(|(x, y)| *x += y);
                    cur_hap.n += 1;
                    // Partials are temporary, so we don't need to do this
                    // cur_hap.samples_idx |= k.1.samples_idx | k.0.
                }
                cur_hap.partial = m_len - i;
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
        let coverage_ordering = self.meta.coverage.cmp(&other.meta.coverage);
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

        self.kfeat
            .iter()
            .zip(&other.kfeat)
            .find_map(|(i, j)| {
                if (*i as u64) != (*j as u64) {
                    Some((*i as u64).cmp(&(*j as u64)))
                } else {
                    None
                }
            })
            .unwrap_or(Ordering::Equal)
    }
}

impl PartialEq for Haplotype {
    fn eq(&self, other: &Self) -> bool {
        self.meta.coverage == other.meta.coverage
            && self.size == other.size
            && self.n == other.n
            && self
                .kfeat
                .iter()
                .zip(&other.kfeat)
                .all(|(i, j)| *i as u64 == *j as u64)
    }
}

impl Eq for Haplotype {}

impl Hash for Haplotype {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for &val in &self.kfeat {
            val.to_bits().hash(state);
        }
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
            .field("samp", &self.meta.samples_idx)
            // Exclude kfeat from the debug output
            .finish()
    }
}

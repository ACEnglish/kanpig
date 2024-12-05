use crate::kplib::seq_to_kmer;
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
    pub coverage: u64,
    pub kfeat: Vec<f32>,
    pub parts: Vec<(i64, Vec<f32>)>,
    pub partial: usize,
    pub ps: Option<u16>,
    pub hp: Option<u8>,
}

impl Haplotype {
    pub fn new(
        kfeat: Vec<f32>,
        size: i64,
        n: u64,
        coverage: u64,
        ps: Option<u16>,
        hp: Option<u8>,
    ) -> Self {
        Self {
            size,
            n,
            coverage,
            kfeat: kfeat.clone(),
            parts: vec![(size, kfeat)],
            partial: 0,
            ps,
            hp,
        }
    }

    // Create an empty haplotype
    pub fn blank(kmer: u8, coverage: u64) -> Haplotype {
        let mk = seq_to_kmer(&[], kmer, false, 0);
        Haplotype {
            size: 0,
            n: 0,
            coverage,
            kfeat: mk.clone(),
            parts: vec![],
            partial: 0,
            ps: None,
            hp: None,
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
                let mut cur_hap = Haplotype::blank(kmer, self.coverage);
                for k in j.iter() {
                    cur_hap.size += k.0;
                    cur_hap
                        .kfeat
                        .iter_mut()
                        .zip(k.1.iter())
                        .for_each(|(x, y)| *x += y);
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
        let coverage_ordering = self.coverage.cmp(&other.coverage);
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
        self.coverage == other.coverage
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
            .field("coverage", &self.coverage)
            .field("ps", &self.ps)
            .field("hp", &self.hp)
            // Exclude kfeat from the debug output
            .finish()
    }
}

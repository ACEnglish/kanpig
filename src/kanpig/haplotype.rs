use crate::kanpig::seq_to_kmer;
use itertools::Itertools;
use std::cmp::Ordering;
use std::fmt::{Debug, Formatter, Result};

#[derive(Clone)]
pub struct Haplotype {
    pub size: i64,
    pub n: u64,
    pub coverage: u64,
    pub kfeat: Vec<f32>,
    pub parts: Vec<(i64, Vec<f32>)>,
}

impl Haplotype {
    pub fn new(kfeat: Vec<f32>, size: i64, n: u64, coverage: u64) -> Self {
        Haplotype {
            size: size,
            n,
            coverage,
            kfeat: kfeat.clone(),
            parts: vec![(size, kfeat)],
        }
    }

    // Create an empty haplotype
    pub fn blank(kmer: u8, coverage: u64) -> Haplotype {
        let mk = seq_to_kmer(&[], kmer, false);
        Haplotype {
            size: 0,
            n: 0,
            coverage,
            kfeat: mk.clone(),
            parts: vec![],
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

    pub fn partial_haplotypes(&self) -> Vec<Haplotype> {
        let mut ret = vec![];
        for i in (1..=self.parts.len()).rev() {
            for j in self.parts.iter().combinations(i) {
                let mut cur_hap = Haplotype::blank(4, self.coverage);
                for k in j.iter() {
                    cur_hap.size += k.0;
                    cur_hap
                        .kfeat
                        .iter_mut()
                        .zip(k.1.iter())
                        .for_each(|(x, y)| *x += y);
                }
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

impl Debug for Haplotype {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        f.debug_struct("Haplotype")
            .field("size", &self.size)
            .field("n", &self.n)
            .field("coverage", &self.coverage)
            // Exclude kfeat from the debug output
            .finish()
    }
}

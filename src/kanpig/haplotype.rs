use crate::kanpig::seq_to_kmer;

#[derive(Clone)]
pub struct Haplotype {
    pub size: i64,
    pub n: u64,
    pub coverage: u64,
    pub gq: Option<Vec<f64>>,
    pub kfeat: Vec<f32>,
}

impl Haplotype {
    pub fn new(kfeat: Vec<f32>, size: i64, n: u64, coverage: u64) -> Self {
        Haplotype {
            size,
            n,
            coverage,
            gq: None,
            kfeat,
        }
    }

    // Create an empty haplotype
    pub fn blank(kmer: u8, coverage: u64) -> Haplotype {
        Haplotype {
            size: 0,
            n: 0,
            coverage,
            gq: None,
            kfeat: seq_to_kmer(&[], kmer, false),
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
    }
}

impl PartialOrd for Haplotype {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Haplotype {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // First, compare by coverage. More coverage is better
        let coverage_ordering = self.coverage.partial_cmp(&other.coverage).unwrap();
        if coverage_ordering != std::cmp::Ordering::Equal {
            return coverage_ordering;
        }

        // If coverage is equal, compare by number of variants
        // We prefer fewer variants, thus the reverse
        let changes_ordering = self.n.partial_cmp(&other.n).unwrap();
        if changes_ordering != std::cmp::Ordering::Equal {
            return changes_ordering.reverse(); 
        }
        // sort by size - This makes a preference for keeping smaller SVs
        self.size.cmp(&other.size)
    }
}

impl PartialEq for Haplotype {
    fn eq(&self, other: &Self) -> bool {
        let first = self.coverage == other.coverage && self.size == other.size && self.n == other.n;
        if !first {
            return false;
        }
        for (i, j) in self.kfeat.iter().zip(other.kfeat.iter()) {
            if *i as u64 != *j as u64 {
                return false;
            }
        }
        true
    }
}

impl Eq for Haplotype {}


impl std::fmt::Debug for Haplotype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Haplotype")
            .field("size", &self.size)
            .field("n", &self.n)
            .field("coverage", &self.coverage)
            .field("gq", &self.gq)
            // Exclude kfeat from the debug output
            .finish()
    }
}


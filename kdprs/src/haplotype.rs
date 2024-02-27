use crate::kmer::seq_to_kmer;

#[derive(Clone, Debug)]
pub struct Haplotype {
    pub kfeat: Vec<f32>,
    pub size: i64,
    pub n: u64,
    pub coverage: u64,
}

impl Haplotype {
    pub fn new(kfeat: Vec<f32>, size: i64, n: u64, coverage: u64) -> Self {
        Haplotype {
            kfeat,
            size,
            n,
            coverage,
        }
    }

    // Create an empty haplotype
    pub fn blank(kmer: u8, coverage: u64) -> Haplotype {
        Haplotype {
            kfeat: seq_to_kmer(&[], kmer),
            size: 0,
            n: 0,
            coverage,
        }
    }

    // Add another variant to a Haplotype
    pub fn add(&mut self, other: &Haplotype) {
        if !self.kfeat.len() == other.kfeat.len() {
            panic!("Cannot add haplotypes of different kmer size");
        }
        let _ = self
            .kfeat
            .iter_mut()
            .zip(other.kfeat.iter())
            .for_each(|(x, y)| *x += y);
        self.size += other.size;
        self.n += 1;
    }
}

impl PartialOrd for Haplotype {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        // First, compare by coverage
        let coverage_ordering = self.coverage.partial_cmp(&other.coverage)?;
        if coverage_ordering != std::cmp::Ordering::Equal {
            return Some(coverage_ordering);
        }

        // If coverage is equal, compare by number of variants
        Some(self.n.cmp(&other.n))
    }
}

impl Ord for Haplotype {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
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

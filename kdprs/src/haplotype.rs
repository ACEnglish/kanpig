use crate::kmer::seq_to_kmer;

#[derive(Debug, Clone)]
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
            .zip(&other.kfeat)
            .map(|(x, y)| *x += y);
        self.size += other.size;
        self.n += 1;
    }
}

// Can I do a ::new but also call it directly via Haplotype { a, b, c,.. }

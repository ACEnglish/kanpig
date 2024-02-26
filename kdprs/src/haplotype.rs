use crate::kmer::seq_to_kmer;

pub struct Haplotype {
    kfeat: Vec<f32>,
    size: i64,
    n: u64,
    coverage: u64,
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
}

// Can I do a ::new but also call it directly via Haplotype { a, b, c,.. }

/// A pileup variant that's hashable / comparable
use crate::kplib::Svtype;
use rust_htslib::bam::pileup::{Alignment, Indel};
use std::hash::{Hash, Hasher};

#[derive(Clone)]
pub struct PileupVariant {
    pub position: u64,
    pub indel: Svtype,
    pub size: i64,
    pub sequence: Option<Vec<u8>>,
}

impl PileupVariant {
    pub fn new(alignment: &Alignment, position: u64) -> Self {
        let (indel, size, sequence) = match alignment.indel() {
            Indel::Del(size) => (Svtype::Del, -(size as i64), None),
            Indel::Ins(size) => {
                let qpos = alignment.qpos().unwrap();
                let seq =
                    alignment.record().seq().as_bytes()[qpos..(qpos + size as usize)].to_vec();
                (Svtype::Ins, size as i64, Some(seq))
            }
            _ => panic!("Unexpected Indel type"),
        };

        PileupVariant {
            position,
            indel,
            size,
            sequence,
        }
    }
}

// Implement PartialEq trait for custom equality comparison
impl PartialEq for PileupVariant {
    fn eq(&self, other: &Self) -> bool {
        // Compare position, size, and indel types
        if self.position != other.position || self.size != other.size || self.indel != other.indel {
            return false;
        }

        // If the indel is a deletion, only compare positions and sizes
        if self.indel == Svtype::Del {
            return true;
        }

        // Otherwise, compare sequences
        self.sequence == other.sequence
    }
}

// Implement Eq trait for reflexive equality
impl Eq for PileupVariant {}

impl Hash for PileupVariant {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.position.hash(state);
        self.size.hash(state);
        self.indel.hash(state);
        if self.indel == Svtype::Ins {
            self.sequence.as_ref().unwrap().hash(state);
        }
    }
}

impl std::fmt::Debug for PileupVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let seq = match &self.sequence {
            Some(seq) => String::from_utf8_lossy(seq),
            None => String::from("").into(),
        };
        f.debug_struct("PileupVariant")
            .field("position", &self.position)
            .field("size", &self.size)
            .field("indel", &self.indel)
            .field("sequence", &seq)
            // Exclude kfeat from the debug output
            .finish()
    }
}

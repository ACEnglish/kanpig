/*
 * A pileup variant that's hashable / equatable
 */
use crate::kanplug::Svtype;
use rust_htslib::bam::pileup::{Alignment, Indel};
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone)]
pub struct PileupVariant {
    pub position: u64,
    pub size: i64,
    pub sequence: Option<Vec<u8>>,
    pub indel: Svtype,
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
            _ => panic!("this can't happen"),
        };

        PileupVariant {
            position,
            size,
            sequence,
            indel,
        }
    }
    // pub fn to_hap(&self) -> Haplotype going to need reference
}

// Implement PartialEq trait for custom equality comparison
impl PartialEq for PileupVariant {
    fn eq(&self, other: &Self) -> bool {
        let first =
            self.position == other.position && self.size == other.size && self.indel == other.indel;
        if !first {
            return false;
        }
        // Don't compare DEL sequences
        if self.indel == Svtype::Del {
            true
        } else {
            self.sequence == other.sequence
        }
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

/// A pileup variant that's hashable / comparable
use crate::kplib::Svtype;
use crate::kplib::{seq_to_kmer, Haplotype, KDParams};
use indexmap::{IndexMap, IndexSet};
use rust_htslib::faidx;
use std::hash::{Hash, Hasher};

#[derive(Clone)]
pub struct PileupVariant {
    pub position: u64,
    pub end: u64,
    pub indel: Svtype,
    pub size: i64,
    pub sequence: Option<Vec<u8>>,
}

impl PileupVariant {
    pub fn new(
        position: u64,
        end: u64,
        indel: Svtype,
        size: i64,
        sequence: Option<Vec<u8>>,
    ) -> Self {
        PileupVariant {
            position,
            end,
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

pub type ReadsMap = IndexMap<Vec<u8>, Vec<usize>>;
pub type PileupSet = IndexSet<PileupVariant>;

pub fn pileups_to_haps(
    chrom: &String,
    reads: ReadsMap,
    mut plups: PileupSet,
    coverage: u64,
    reference: &faidx::Reader,
    params: &KDParams,
) -> (Vec<Haplotype>, u64) {
    let mut hap_parts = Vec::<Haplotype>::with_capacity(plups.len());

    while let Some(mut p) = plups.pop() {
        // Need to fill in deleted sequence
        let sequence = if p.indel == Svtype::Del {
            let d_start = p.position;
            let d_end = d_start + p.size.unsigned_abs();
            reference
                .fetch_seq(chrom, d_start as usize, d_end as usize)
                .unwrap()
                .to_vec()
        } else {
            p.sequence
                .take()
                .expect("Insertions should already have a sequence")
        };

        let n_hap = Haplotype::new(
            seq_to_kmer(
                &sequence,
                params.kmer,
                p.indel == Svtype::Del,
                params.maxhom,
            ),
            p.size,
            1,
            1,
        );
        hap_parts.push(n_hap);
    }

    // Deduplicate reads by pileup combination
    let mut unique_reads: IndexMap<&Vec<usize>, u64> = IndexMap::new();
    for m_plups in reads.values() {
        *unique_reads.entry(m_plups).or_insert(0) += 1;
    }

    // Turn variants into haplotypes
    let mut ret = Vec::<Haplotype>::new();
    for (read_pileups, coverage) in unique_reads {
        let mut cur_hap = Haplotype::blank(params.kmer, 1);
        for p in read_pileups {
            cur_hap.add(&hap_parts[hap_parts.len() - *p - 1]);
        }
        cur_hap.coverage = coverage;
        ret.push(cur_hap);
    }

    ret.sort_by(|a, b| b.cmp(a));
    trace!("{:?}", ret);
    (ret, coverage)
}

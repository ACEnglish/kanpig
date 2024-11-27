/// A pileup variant that's hashable / comparable
use crate::kplib::Svtype;
use rust_htslib::{bam::ext::BamRecordExtensions, bam::Record};
use std::hash::{Hash, Hasher};

pub struct ReadPileup {
    pub chrom: i32,
    pub start: u64,
    pub end: u64,
    pub pileups: Vec<PileupVariant>,
}

impl ReadPileup {
    pub fn new(record: Record, sizemin: u32, sizemax: u32) -> Self {
        let chrom = record.tid();
        let start = record.reference_start();
        let end = record.reference_end();

        let mut pileups = Vec::<PileupVariant>::new();
        let mut read_offset = 0;
        let mut align_offset = start as usize - 1;

        for cigar in record.cigar().iter() {
            match cigar.char() {
                'I' if sizemin <= cigar.len() && cigar.len() <= sizemax => {
                    let sequence = record.seq().as_bytes()
                        [read_offset..read_offset + cigar.len() as usize]
                        .to_vec();

                    pileups.push(PileupVariant::new(
                        align_offset as u64,
                        (align_offset + 1) as u64,
                        Svtype::Ins,
                        cigar.len() as i64,
                        Some(sequence),
                    ));

                    read_offset += cigar.len() as usize;
                }
                'D' if sizemin <= cigar.len() && cigar.len() <= sizemax => {
                    pileups.push(PileupVariant::new(
                        align_offset as u64,
                        (align_offset + (cigar.len() as usize)) as u64,
                        Svtype::Del,
                        -(cigar.len() as i64),
                        None,
                    ));
                    align_offset += cigar.len() as usize;
                }
                'M' | 'X' | '=' => {
                    read_offset += cigar.len() as usize;
                    align_offset += cigar.len() as usize;
                }
                'I' | 'S' => {
                    read_offset += cigar.len() as usize;
                }
                'D' => {
                    align_offset += cigar.len() as usize;
                }
                // Handle 'H' and 'P' explicitly to ignore them
                'H' | 'P' => {}
                _ => {
                    error!("What is this code?: {}", cigar.char());
                }
            }
        }

        ReadPileup {
            chrom,
            start: start as u64,
            end: end as u64,
            pileups,
        }
    }

    pub fn decode(line: &[u8], sizemin: u64, sizemax: u64) -> Option<Self> {
        let line_str = std::str::from_utf8(line).ok()?;
        let mut fields = line_str.split('\t');

        let _chrom = fields.next()?.to_string();
        let start = fields.next()?.parse().ok()?;
        let end = fields.next()?.parse().ok()?;
        let pileups_str = fields.next()?;

        let pileups = match pileups_str {
            "." => Vec::new(),
            _ => pileups_str
                .split(',')
                .filter_map(|entry| PileupVariant::decode(entry, start))
                .filter(|variant| variant.size.unsigned_abs() >= sizemin && variant.size.unsigned_abs() <= sizemax)
                .collect(),
        };
        // I use chrom 0 for the decode because new puts in tid
        Some(ReadPileup {
            chrom: 0,
            start,
            end,
            pileups,
        })
    }

    pub fn to_string(&self, chrom: &str) -> String {
        let pstr: String = if self.pileups.is_empty() {
            ".".to_string()
        } else {
            self.pileups
                .iter()
                .map(|x| x.encode(self.start))
                .collect::<Vec<_>>()
                .join(",")
        };
        format!("{}\t{}\t{}\t{}", chrom, self.start, self.end, pstr)
    }
}

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

    /// Parses a single entry into a `PileupVariant`.
    pub fn decode(entry: &str, start: u64) -> Option<Self> {
        let mut parts = entry.split(':');
        let offset = parts.next()?.parse::<u64>().ok()?;
        let m_pos = start + offset;
        let value = parts.next()?;

        let (end, svtype, size, seq) = match value.parse::<u64>() {
            Ok(size) => (m_pos + size, Svtype::Del, -(size as i64), None),
            Err(_) => (
                m_pos + 1,
                Svtype::Ins,
                value.len() as i64,
                Some(value.as_bytes().to_vec()),
            ),
        };

        Some(PileupVariant::new(m_pos, end, svtype, size, seq))
    }

    pub fn encode(&self, offset: u64) -> String {
        match self.indel {
            Svtype::Del => format!("{}:{}", self.position - offset, self.size.abs()),
            Svtype::Ins => format!(
                "{}:{}",
                self.position - offset,
                self.sequence
                    .clone()
                    .and_then(|seq| String::from_utf8(seq).ok())
                    .unwrap()
            ),
            _ => panic!("Unencodeable PileupVariant"),
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

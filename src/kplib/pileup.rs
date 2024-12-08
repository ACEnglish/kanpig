/// A pileup variant that's hashable / comparable
use crate::kplib::Svtype;
use rust_htslib::{bam::ext::BamRecordExtensions, bam::record::Aux, bam::Record};
use std::hash::{Hash, Hasher};

#[derive(Debug)]
pub struct ReadPileup {
    pub chrom: i32,
    pub start: u64,
    pub end: u64,
    pub pileups: Vec<PileupVariant>,
    pub ps: Option<u16>,
    pub hp: Option<u8>,
}

/// A struct representing a read and its pileups
impl ReadPileup {
    /// Creates a new `ReadPileup` from an alignment record, filtering variants based on size constraints.
    ///
    /// # Parameters
    /// - `record`: The `Record` from which the pileup is constructed.
    /// - `sizemin`: The minimum size of variants to include in the pileup.
    /// - `sizemax`: The maximum size of variants to include in the pileup.
    ///
    /// # Returns
    /// - A `ReadPileup` instance containing the extracted variants that satisfy the size constraints.
    ///
    /// # Function Details
    /// - Iterates through the record's CIGAR operations and processes the following:
    ///     - **Insertions (I)**: Captures the inserted sequence if its length is within `[sizemin, sizemax]`.
    ///     - **Deletions (D)**: Captures the deletion if its length is within `[sizemin, sizemax]`.
    ///     - **Matches (M, X, =) and Soft-clips (S)**: Advances offsets without creating variants.
    ///     - **Hard-clips (H) and Pads (P)**: Ignores these operations.
    /// - Logs an error for any unexpected CIGAR operation.
    ///
    /// # Example
    /// ```rust
    /// let record = ...; // A valid BAM record
    /// let pileup = ReadPileup::new(record, 10, 100);
    /// println!("{:?}", pileup);
    /// ```
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

        let ps = match record.aux(b"PS") {
            Ok(Aux::U16(value)) => Some(value),
            _ => None,
        };

        let hp = match record.aux(b"HP") {
            Ok(Aux::U8(value)) => Some(value),
            _ => None,
        };

        Self {
            chrom,
            start: start as u64,
            end: end as u64,
            pileups,
            ps,
            hp,
        }
    }

    /// Decodes a `ReadPileup` from a tab-delimited string, applying size constraints to the included variants.
    ///
    /// # Parameters
    /// - `line`: A byte slice containing the tab-delimited representation of a pileup.
    /// - `sizemin`: The minimum size of variants to include in the pileup.
    /// - `sizemax`: The maximum size of variants to include in the pileup.
    ///
    /// # Returns
    /// - `Option<Self>`: Returns `Some(ReadPileup)` if decoding is successful, or `None` if the input is invalid or contains no valid variants.
    ///
    /// # Function Details
    /// - Splits the input line into fields:
    ///     - `chrom`: The chromosome name (ignored during decoding, set to 0).
    ///     - `start` and `end`: The start and end positions of the pileup.
    ///     - `pileups`: A comma-separated string of variants.
    /// - Filters decoded variants based on size constraints.
    ///
    /// # TODO: Since parsing from an alignment record uses tid (i32) for chrom, the chromosome is
    /// never parsed from a string.
    ///
    /// # Example
    /// ```rust
    /// let line = b"chr1\t1000\t1010\t.";
    /// let pileup = ReadPileup::decode(line, 10, 100);
    /// println!("{:?}", pileup);
    /// ```
    pub fn decode(line: &[u8], sizemin: u32, sizemax: u32) -> Option<Self> {
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
                .filter(|variant| {
                    variant.size.unsigned_abs() >= sizemin as u64
                        && variant.size.unsigned_abs() <= sizemax as u64
                })
                .collect(),
        };

        let ps = fields.next()?;
        let ps = match ps {
            "." => None,
            _ => Some(ps.parse().ok()?),
        };

        let hp = fields.next()?;
        let hp = match hp {
            "." => None,
            _ => Some(hp.parse().ok()?),
        };

        // I use chrom 0 for the decode because new puts in tid
        Some(ReadPileup {
            chrom: 0,
            start,
            end,
            pileups,
            ps,
            hp,
        })
    }

    /// Encodes a `ReadPileup` into a tab-delimited string representation.
    ///
    /// # Parameters
    /// - `chrom`: The chromosome name as a `&str`.
    ///
    /// # Returns
    /// - `String`: A tab-delimited string containing:
    ///     - Chromosome name
    ///     - Start position
    ///     - End position
    ///     - A comma-separated list of encoded PileupVariants, or `.` if there are no variants.
    ///
    /// # Example
    /// ```rust
    /// let pileup = ReadPileup { chrom: 1, start: 1000, end: 1010, pileups: Vec::new() };
    /// let pileup_str = pileup.to_string("chr1");
    /// println!("{}", pileup_str); // Output: "chr1\t1000\t1010\t."
    /// ```
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

        let ps: String = match self.ps {
            Some(p) => p.to_string(),
            None => ".".to_string(),
        };

        let hp: String = match self.hp {
            Some(h) => h.to_string(),
            None => ".".to_string(),
        };

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            chrom, self.start, self.end, pstr, ps, hp
        )
    }
}

pub struct PileupVariant {
    pub position: u64,
    pub end: u64,
    pub indel: Svtype,
    pub size: i64,
    pub sequence: Option<Vec<u8>>,
}

/// Provides information for an individual deletion or insertion with
/// methods to create, decode, and encode structural variant information.
impl PileupVariant {
    /// Creates a new `PileupVariant`.
    ///
    /// # Parameters
    /// - `position`: The start position of the variant.
    /// - `end`: The end position of the variant.
    /// - `indel`: The type of the structural variant (`Svtype::Del` or `Svtype::Ins`).
    /// - `size`: The size of the variant (positive for insertions, negative for deletions).
    /// - `sequence`: The sequence of the variant (only applicable for insertions).
    ///
    /// # Returns
    /// A new `PileupVariant` instance with the specified properties.
    ///
    /// # Example
    /// ```rust
    /// let variant = PileupVariant::new(1000, 1001, Svtype::Ins, 50, Some(vec![65, 67, 71, 84]));
    /// ```
    pub fn new(
        position: u64,
        end: u64,
        indel: Svtype,
        size: i64,
        sequence: Option<Vec<u8>>,
    ) -> Self {
        Self {
            position,
            end,
            indel,
            size,
            sequence,
        }
    }

    /// Decodes a string entry into a `PileupVariant`.
    ///
    /// # Parameters
    /// - `entry`: A string slice representing a variant entry (e.g., `offset:size` for deletions or `offset:sequence` for insertions).
    /// - `start`: The start position of the reference region to calculate the absolute position.
    ///
    /// # Returns
    /// - `Some(PileupVariant)` if the entry is valid.
    /// - `None` if the entry cannot be parsed.
    ///
    /// # Parsing Logic
    /// - For deletions (`offset:size`):
    ///   - Parses `offset` and `size` to calculate the variant position and size.
    ///   - Sets `sequence` to `None`.
    /// - For insertions (`offset:sequence`):
    ///   - Parses `offset` and derives the sequence from the remaining string.
    ///   - Calculates the size based on the sequence length.
    ///
    /// # Example
    /// ```rust
    /// let variant = PileupVariant::decode("10:ACGT", 1000).unwrap();
    /// assert_eq!(variant.position, 1010);
    /// assert_eq!(variant.indel, Svtype::Ins);
    /// assert_eq!(variant.size, 4);
    /// ```
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

    /// Encodes a `PileupVariant` into a string representation.
    ///
    /// # Parameters
    /// - `offset`: The reference start position to calculate relative positions.
    ///
    /// # Returns
    /// A string representation of the variant:
    /// - For deletions: `offset:size` (e.g., `10:50`).
    /// - For insertions: `offset:sequence` (e.g., `10:ACGT`).
    ///
    /// # Panics
    /// - If the variant type is neither `Svtype::Del` nor `Svtype::Ins`.
    ///
    /// # Example
    /// ```rust
    /// let variant = PileupVariant::new(1010, 1011, Svtype::Ins, 4, Some(vec![65, 67, 71, 84]));
    /// assert_eq!(variant.encode(1000), "10:ACGT");
    /// ```
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

use crate::kplib::{KDParams, PileupSet, PileupVariant, ReadParser, ReadsMap};
use rust_htslib::tbx::{self, Read};
use std::path::PathBuf;

pub struct PlupParser {
    tbx: tbx::Reader,
    params: KDParams,
}

impl PlupParser {
    /// Creates a new `PlupReader` for a given file path.
    pub fn new(file_path: PathBuf, params: KDParams) -> Self {
        let tbx = tbx::Reader::from_path(&file_path).expect("Failed to open TBX file");
        Self { tbx, params }
    }

    /// Parses a single line from the file.
    fn parse_line(line: &[u8]) -> Option<(String, u64, u64, Vec<PileupVariant>)> {
        let line_str = std::str::from_utf8(line).ok()?;
        let mut fields = line_str.split('\t');

        let chrom = fields.next()?.to_string();
        let start = fields.next()?.parse().ok()?;
        let end = fields.next()?.parse().ok()?;
        let pileups_str = fields.next()?;

        let pileups = match pileups_str {
            "." => Vec::new(),
            _ => pileups_str
                .split(',')
                .filter_map(|entry| PileupVariant::decode(entry, start))
                .collect(),
        };

        Some((chrom, start, end, pileups))
    }
}

impl ReadParser for PlupParser {
    /// Fetch and parse pileups within a specified genomic interval.
    fn find_pileups(&mut self, chrom: &str, start: u64, end: u64) -> (ReadsMap, PileupSet, u64) {
        let window_start = start.saturating_sub(self.params.chunksize);
        let window_end = end + self.params.chunksize;

        let tid = self
            .tbx
            .tid(chrom)
            .unwrap_or_else(|_| panic!("Could not resolve '{}' to contig ID", chrom));
        self.tbx
            .fetch(tid, window_start, window_end)
            .expect("Could not fetch region from TBX");

        let mut reads = ReadsMap::new();
        let mut p_variants = PileupSet::new();
        let mut total_coverage = 0;

        for (qname, line) in self.tbx.records().filter_map(Result::ok).enumerate() {
            if let Some((_, start, end, variants)) = Self::parse_line(&line) {
                if start <= window_start && end >= window_end {
                    total_coverage += 1;
                    for m_var in variants {
                        let size = m_var.size.unsigned_abs();
                        if m_var.position > window_start
                            && m_var.position < window_end
                            && size >= self.params.sizemin
                            && size <= self.params.sizemax
                        {
                            trace!("{:?}", m_var);
                            let (p_idx, _) = p_variants.insert_full(m_var);
                            reads.entry(qname).or_default().push(p_idx);
                        }
                    }
                }
            }
        }

        (reads, p_variants, total_coverage)
    }
}

use crate::kplib::{KDParams, PileupSet, PileupVariant, ReadParser, ReadsMap, Svtype};
use rust_htslib::tbx::{self, Read};
use std::path::PathBuf;

pub struct PlupParser {
    tbx: tbx::Reader,
    params: KDParams,
}

impl PlupParser {
    /// Creates a new `PlupReader` for a given file path.
    pub fn new(file_path: PathBuf, params: KDParams) -> Self {
        let tbx = tbx::Reader::from_path(&file_path).unwrap();
        Self { tbx, params }
    }

    /// Parses a single line from the file.
    fn parse_line(line: &[u8]) -> Option<(String, u64, u64, Vec<PileupVariant>)> {
        let line_str = std::str::from_utf8(line).ok()?;
        let mut fields = line_str.split('\t');

        let chrom = fields.next()?.to_string();
        let start: u64 = fields.next()?.parse().ok()?;
        let end: u64 = fields.next()?.parse().ok()?;
        let pileups_str = fields.next()?;

        let pileups = pileups_str
            .split(',')
            .filter_map(|entry| {
                let mut parts = entry.split(':');
                let offset: u64 = parts.next()?.parse().ok()?;
                let m_pos = start + offset;
                let value = parts.next()?;
                let (end, svtype, size, seq) = if let Ok(size) = value.parse::<u64>() {
                    // It's an integer, so should be Del
                    (m_pos + size, Svtype::Del, -(size as i64), None)
                } else {
                    (
                        m_pos + 1,
                        Svtype::Ins,
                        value.len() as i64,
                        Some(value.as_bytes().to_vec()),
                    )
                };

                Some(PileupVariant::new(m_pos - 1, end, svtype, size, seq))
            })
            .collect();

        Some((chrom, start, end, pileups))
    }
}

impl ReadParser for PlupParser {
    /// Fetch and parse pileups within a specified genomic interval.
    fn find_pileups(&mut self, chrom: &str, start: u64, end: u64) -> (ReadsMap, PileupSet, u64) {
        let window_start = if start < self.params.chunksize {
            0
        } else {
            start - self.params.chunksize
        };
        let window_end = end + self.params.chunksize;

        let tid = match self.tbx.tid(chrom) {
            Ok(tid) => tid,
            Err(_) => panic!("Could not resolve '{}' to contig ID", chrom),
        };
        // Set region to fetch.
        self.tbx
            .fetch(tid, window_start, window_end)
            .expect("Could not seek to {}:{}-{}");

        // track the changes made by each read
        let mut reads = ReadsMap::new();
        // consolidate common variants
        let mut p_variants = PileupSet::new();
        let mut tot_cov: u64 = 0;

        let mut qname = 0;
        for line in self.tbx.records().filter_map(Result::ok) {
            if let Some(parsed) = Self::parse_line(&line) {
                if !(parsed.1 < window_start && parsed.2 > window_end) {
                    continue;
                }
                tot_cov += 1;
                for m_var in parsed.3 {
                    let lsize = m_var.size.unsigned_abs();
                    if m_var.position > window_start
                        && m_var.position < window_end
                        && lsize >= self.params.sizemin
                        && lsize <= self.params.sizemax
                    {
                        let (p_idx, _) = p_variants.insert_full(m_var);
                        reads
                            .entry(qname.to_string().into_bytes())
                            .or_default()
                            .push(p_idx);
                    }
                }
            }
            qname += 1;
        }
        let coverage = (tot_cov / (window_end - window_start)).max(reads.len() as u64);
        (reads, p_variants, coverage)
    }
}

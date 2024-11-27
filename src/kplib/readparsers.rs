use crate::kplib::{KDParams, PileupSet, ReadPileup, ReadsMap};
use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::{IndexedReader, Read as BamRead},
    tbx::{self, Read as TbxRead},
};
use std::path::PathBuf;

pub trait ReadParser {
    fn find_pileups(&mut self, chrom: &str, start: u64, end: u64) -> (ReadsMap, PileupSet, u64);
}

pub struct BamParser {
    bam: IndexedReader,
    params: KDParams,
}

impl BamParser {
    pub fn new(bam_name: PathBuf, ref_name: PathBuf, params: KDParams) -> Self {
        let mut bam = IndexedReader::from_path(bam_name).unwrap();
        let _ = bam.set_reference(ref_name.clone());
        BamParser { bam, params }
    }
}

impl ReadParser for BamParser {
    /// Returns all unique haplotypes over a region
    fn find_pileups(&mut self, chrom: &str, start: u64, end: u64) -> (ReadsMap, PileupSet, u64) {
        // We pileup a little outside the region for variants
        let window_start = start.saturating_sub(self.params.chunksize);
        let window_end = end + self.params.chunksize;

        if let Err(e) = self.bam.fetch((&chrom, window_start, window_end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        // track the changes made by each read
        let mut reads = ReadsMap::new();
        // consolidate common variants
        let mut p_variants = PileupSet::new();
        let mut coverage = 0;

        for (qname, record) in self
            .bam
            .records()
            .filter_map(|r| {
                r.ok().filter(|rec| {
                    !rec.seq().is_empty()
                        && rec.mapq() >= self.params.mapq
                        && (rec.flags() & self.params.mapflag) == 0
                        && (rec.reference_start() as u64) < window_start
                        && (rec.reference_end() as u64) > window_end
                })
            })
            .enumerate()
        {
            coverage += 1;
            let mut m_plups = ReadPileup::new(
                record,
                self.params.sizemin as u32,
                self.params.sizemax as u32,
            );
            for m_var in m_plups.pileups.drain(..) {
                if m_var.position >= window_start && m_var.position <= window_end {
                    trace!("{:?}", m_var);
                    let (p_idx, _) = p_variants.insert_full(m_var);
                    reads.entry(qname).or_default().push(p_idx);
                }
            }
        }
        (reads, p_variants, coverage)
    }
}

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
            if let Some(read) = ReadPileup::decode(&line) {
                if read.start <= window_start && read.end >= window_end {
                    total_coverage += 1;
                    for m_var in read.pileups {
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

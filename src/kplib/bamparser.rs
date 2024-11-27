use crate::kplib::{KDParams, PileupSet, ReadParser, ReadPileup, ReadsMap};

use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::{IndexedReader, Read},
};
use std::path::PathBuf;

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

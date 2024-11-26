use crate::kplib::{KDParams, PileupSet, PileupVariant, ReadParser, ReadsMap, Svtype};

use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::pileup::Indel,
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
        let window_start = if start < self.params.chunksize {
            0
        } else {
            start - self.params.chunksize
        };
        let window_end = end + self.params.chunksize;

        let mut tot_cov: u64 = 0;
        if let Err(e) = self.bam.fetch((&chrom, window_start, window_end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        // track the changes made by each read
        let mut reads = ReadsMap::new();
        // consolidate common variants
        let mut p_variants = PileupSet::new();

        for pileup in self.bam.pileup().filter_map(Result::ok) {
            let m_pos: u64 = pileup.pos().into();

            // We got to truncate the pileup
            if m_pos < window_start {
                continue;
            }
            // Do this separately so we don't waste time downstream
            if window_end < m_pos {
                break;
            }

            // Only count coverage of reads used?
            //tot_cov += pileup.depth() as u64;
            for alignment in pileup.alignments() {
                // Skip records without sequence, below min mapq, matching the flag,
                // or partially into our window. Will revisit partial when I can turn a Haplotype into a
                // single-path graph
                if alignment.record().seq().is_empty()
                    || alignment.record().mapq() < self.params.mapq
                    || (alignment.record().flags() & self.params.mapflag) != 0
                    || (self.params.spanoff
                        && !((alignment.record().reference_start() as u64) < window_start
                            && (alignment.record().reference_end() as u64) > window_end))
                {
                    continue;
                }

                // Only count reads we're using - do this before we worry about the
                // is_del/is_refskip
                tot_cov += 1;

                // None if either is_del or is_refskip, we we don't need it
                // We do this separately from above so we can count 'deleted'  as covered
                if alignment.qpos().is_none() {
                    continue;
                }

                let (end, size, svtype, seq) = match alignment.indel() {
                    Indel::Del(size) if size as u64 >= self.params.sizemin => {
                        (m_pos + size as u64, -(size as i64), Svtype::Del, None)
                    }
                    Indel::Ins(size) if size as u64 >= self.params.sizemin => {
                        let qpos = alignment.qpos().unwrap() + 1;
                        let seq = alignment.record().seq().as_bytes()[qpos..(qpos + size as usize)]
                            .to_vec();
                        (m_pos + 1_u64, size as i64, Svtype::Ins, Some(seq))
                    }
                    _ => continue,
                };
                let m_var = PileupVariant::new(m_pos, end, svtype, size, seq);
                trace!("{:?}", m_var);

                let (p_idx, _) = p_variants.insert_full(m_var);
                let qname = alignment.record().qname().to_owned();
                reads.entry(qname).or_default().push(p_idx);
            }
        }
        // Either all the reads used or the mean coverage over the window
        let coverage = (tot_cov / (window_end - window_start)).max(reads.len() as u64);

        (reads, p_variants, coverage)
    }
}

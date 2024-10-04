use crate::kplib::{seq_to_kmer, Haplotype, KDParams, PileupVariant, Svtype};
use indexmap::{IndexMap, IndexSet};

use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::pileup::Indel,
    bam::{IndexedReader, Read},
    faidx,
};
use std::path::PathBuf;

pub struct BamParser {
    bam: IndexedReader,
    reference: faidx::Reader,
    params: KDParams,
}

impl BamParser {
    pub fn new(bam_name: PathBuf, ref_name: PathBuf, params: KDParams) -> Self {
        let mut bam = IndexedReader::from_path(bam_name).unwrap();
        let _ = bam.set_reference(ref_name.clone());
        let reference = faidx::Reader::from_path(ref_name).unwrap();
        BamParser {
            bam,
            reference,
            params,
        }
    }

    /// Returns all unique haplotypes over a region
    pub fn find_haps(&mut self, chrom: &String, start: u64, end: u64) -> (Vec<Haplotype>, u64) {
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
        let mut reads = IndexMap::<Vec<u8>, Vec<usize>>::new();
        // consolidate common variants
        let mut p_variants = IndexSet::<PileupVariant>::new();

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

                let m_var = match alignment.indel() {
                    Indel::Ins(size) | Indel::Del(size) if size as u64 >= self.params.sizemin => {
                        PileupVariant::new(&alignment, m_pos)
                    }
                    _ => continue,
                };
                debug!("{:?}", m_var);

                let (p_idx, _) = p_variants.insert_full(m_var);
                let qname = alignment.record().qname().to_owned();
                reads.entry(qname).or_default().push(p_idx);
            }
        }

        // Either all the reads used or the mean coverage over the window
        let coverage = (tot_cov / (window_end - window_start)).max(reads.len() as u64);

        let mut hap_parts = Vec::<Haplotype>::with_capacity(p_variants.len());

        //p_variants.reverse();
        while let Some(mut p) = p_variants.pop() {
            // Need to fill in deleted sequence
            let sequence = if p.indel == Svtype::Del {
                let d_start = p.position;
                let d_end = d_start + p.size.unsigned_abs();
                self.reference
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
                    self.params.kmer,
                    p.indel == Svtype::Del,
                    self.params.maxhom,
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
            let mut cur_hap = Haplotype::blank(self.params.kmer, 1);
            for p in read_pileups {
                cur_hap.add(&hap_parts[hap_parts.len() - *p - 1]);
                // When using p_variants.reverse above reversed in the while pop
                //cur_hap.add(&hap_parts[*p]);
            }
            cur_hap.coverage = coverage;
            //debug!("{:?}", cur_hap);
            ret.push(cur_hap);
        }

        ret.sort_by(|a, b| b.cmp(a));
        debug!("{:?}", ret);
        (ret, coverage)
    }
}

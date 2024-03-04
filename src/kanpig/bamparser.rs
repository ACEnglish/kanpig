use crate::kanpig::{seq_to_kmer, Haplotype, KDParams, PileupVariant, Svtype};
use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::pileup::Indel,
    bam::{IndexedReader, Read},
    faidx,
};
use std::collections::HashMap;
use std::path::PathBuf;

pub struct BamParser {
    bam: IndexedReader,
    reference: faidx::Reader,
    params: KDParams,
}

impl BamParser {
    pub fn new(bam_name: PathBuf, ref_name: PathBuf, params: KDParams) -> Self {
        let bam = IndexedReader::from_path(bam_name).unwrap();
        let reference = faidx::Reader::from_path(ref_name).unwrap();
        BamParser {
            bam,
            reference,
            params,
        }
    }

    // Two possible haplotypes and coverage observed over the window
    pub fn find_haps(&mut self, chrom: &String, start: u64, end: u64) -> (Vec<Haplotype>, u64) {
        // We pileup a little outside the region for variants
        let window_start = start - self.params.chunksize;
        let window_end = end + self.params.chunksize;

        let mut tot_cov: u64 = 0;
        if let Err(e) = self.bam.fetch((
            &chrom,
            start - self.params.chunksize,
            end + self.params.chunksize,
        )) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        // track the changes made by each read
        let mut m_reads = HashMap::<Vec<u8>, Vec<PileupVariant>>::new();
        // consolidate common variants
        let mut p_variants = HashMap::<PileupVariant, u64>::new();
        for pileup in self.bam.pileup() {
            if pileup.is_err() {
                continue;
            }

            let pileup = pileup.unwrap();
            let m_pos: u64 = pileup.pos().into();

            // We got to truncate the pileup
            if m_pos < window_start {
                continue;
            }
            if window_end < m_pos {
                break;
            }

            tot_cov += pileup.depth() as u64;

            for alignment in pileup.alignments() {
                // None if either is_del or is_refskip, we we don't need it
                if alignment.qpos().is_none() {
                    continue;
                }

                // Guard against partial alignments which mess up the kfeat
                // Will revisit when I can turn a Haplotype into a single-path graph
                if !((alignment.record().reference_start() as u64) < start
                    && (alignment.record().reference_end() as u64) > end)
                {
                    continue;
                }

                let m_var = match alignment.indel() {
                    Indel::Ins(size) | Indel::Del(size) if size as u64 >= self.params.sizemin => {
                        PileupVariant::new(&alignment, m_pos)
                    }
                    _ => continue,
                };
                debug!("{:?}", m_var);
                *p_variants.entry(m_var.clone()).or_insert(0) += 1;
                let qname = alignment.record().qname().to_owned();
                m_reads.entry(qname).or_default().push(m_var);
            }
        }

        let m_haps = self.reads_to_haps(m_reads, p_variants, chrom);
        let coverage = tot_cov / (window_end - window_start);
        (m_haps, coverage)
    }

    fn reads_to_haps(
        &self,
        reads: HashMap<Vec<u8>, Vec<PileupVariant>>,
        pileups: HashMap<PileupVariant, u64>,
        chrom: &String,
    ) -> Vec<Haplotype> {
        let mut hap_parts = HashMap::<PileupVariant, Haplotype>::new();

        for (mut p, _) in pileups.into_iter() {
            // Need to fill in deleted sequence
            if p.indel == Svtype::Del {
                let d_start = p.position;
                let d_end = d_start + p.size.unsigned_abs();
                p.sequence = Some(
                    self.reference
                        .fetch_seq(chrom, d_start as usize, d_end as usize)
                        .unwrap()
                        .to_vec(),
                );
            }

            let n_hap = Haplotype::new(
                seq_to_kmer(
                    &p.sequence.clone().unwrap(),
                    self.params.kmer,
                    p.indel == Svtype::Del,
                ),
                p.size,
                1,
                1,
            );
            hap_parts.insert(p, n_hap);
        }

        // Deduplicate reads
        let mut unique_reads: HashMap<&Vec<PileupVariant>, u64> = HashMap::new();
        for m_plups in reads.values() {
            *unique_reads.entry(m_plups).or_insert(0) += 1;
        }

        // Turn variants into haplotypes
        let mut ret = Vec::<Haplotype>::new();
        for (read_pileups, coverage) in unique_reads {
            let mut cur_hap = Haplotype::blank(self.params.kmer, 1);
            for p in read_pileups {
                cur_hap.add(&hap_parts[&p]);
            }
            cur_hap.coverage = coverage;
            debug!("{:?}", cur_hap);
            ret.push(cur_hap);
        }

        // Sort so higher coverage haplotypes are preferred as centroids
        ret.sort_by_key(|i| std::cmp::Reverse(i.coverage));
        ret
    }
}

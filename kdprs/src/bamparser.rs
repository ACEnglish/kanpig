/*
 * Okay, so this thing is going to hold the bam and reference and kdparams
 * Its query is a find_haps method that will take chrom, reg_start, reg_end
 * find_haps will pull the reference sequence (don't always need to...)
 * find_haps will pileup over the regions ±chunksize min_base_quality=0, truncate=true
 * keep track of the total coverage
 * for read in column.pileups:
 *  if query position is none and read.indel not within sizemin/sizemax --
 *  The bam pileup is a little different.. in rust htslib... so going to have to deal with that
 *  but lets just assume we make it all just fine. No optimizations, just a string of Haplotypes
 *
 * We'll also need a hap_deduplicate.. though I think that won't be a thing once we refactor
 * And then we got the noreads or no alts and single read stuff
 * Then there's a read cluster - which isn't too bad with a decent kmeans library
 * just gotta work on the consolidate with best
 */

use crate::cli::KDParams;
use crate::haplotype::Haplotype;
use crate::vcf_traits::Svtype;
use rust_htslib::bam::pileup::{Alignment, Indel};
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::faidx;
use std::collections::HashMap;
use std::path::PathBuf;

pub struct BamParser {
    pub bam_name: PathBuf,
    pub ref_name: PathBuf,
    bam: IndexedReader,
    reference: faidx::Reader,
    params: KDParams,
}

impl BamParser {
    pub fn new(bam_name: PathBuf, ref_name: PathBuf, params: KDParams) -> Self {
        // Open the BAM or CRAM file
        let mut bam = IndexedReader::from_path(&bam_name).unwrap();
        let mut reference = faidx::Reader::from_path(&ref_name).unwrap();
        BamParser {
            bam_name,
            ref_name,
            bam,
            reference,
            params,
        }
    }

    pub fn find_haps(&mut self, chrom: String, start: u64, end: u64) -> (Haplotype, Haplotype) {
        // We pileup a little outside the region for variants
        let window_start = start - self.params.chunksize;
        let window_end = end + self.params.chunksize;

        let mut tot_cov: u64 = 0;
        self.bam.fetch((&chrom, start - self.params.chunksize, end + self.params.chunksize));
        
        let mut m_haps = HashMap::<Vec<u8>, Vec<PileupVariant>>::new();
        for pileup in self.bam.pileup() {
            if !pileup.is_ok() {
                continue;
            }
            
            let pileup = pileup.unwrap();



            let m_pos: u64 = pileup.pos().into();
            // We got to truncate the pileup
            if m_pos < window_start || window_end < m_pos {
                continue;
            }

            tot_cov += pileup.depth() as u64;

            for alignment in pileup.alignments() {
                // None if either is_del or is_refskip, we we don't need it
                if alignment.qpos().is_none() {
                    continue;
                }

                // Guard against partial alignments which mess up the kfeat 
                // Will revisit when I can turn a Haplotype into a single-path graph
                if !((alignment.record().reference_start() as u64) < start && (alignment.record().reference_end() as u64) > end) {
                    continue;
                }

                let m_var = match alignment.indel() {
                    Indel::Ins(size) | Indel::Del(size) if size as u64 >= self.params.sizemin => {
                        PileupVariant::new(&alignment, size, window_start)
                    }
                    _ => continue,
                };

                let qname = alignment.record().qname().to_owned();
                m_haps.entry(qname).or_insert_with(Vec::new).push(m_var);
            }
        }

        let coverage = tot_cov / (window_end - window_start + (2 * self.params.chunksize));
        (Haplotype::blank(3, coverage), Haplotype::blank(3, coverage))
    }

    // fn hap_deduplicate {}
    // fn read_cluster {}
}

struct PileupVariant {
    position: u64,
    size: i64,
    sequence: Option<String>,
    indel: Svtype,
}

impl PileupVariant {
    // https://docs.rs/rust-htslib/latest/rust_htslib/bam/pileup/struct.Alignment.html
    // Need to figure out how to 
    pub fn new(alignment: &Alignment, size: u32, offset: u64) -> Self {
        PileupVariant {
            position: 0,
            size: size as i64,
            sequence: None,
            indel: Svtype::Del,
        }
    }

    // pub fn to_hap(&self) -> Haplotype going to need reference
}

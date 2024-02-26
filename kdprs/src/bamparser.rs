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

use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::bam::pileup::Indel;
use rust_htslib::faidx;
use crate::haplotype::Haplotype;
use crate::cli::KDParams;
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

        // Seek to the region of interest
        self.bam.fetch((&chrom, start, end));

        // Iterate over the pileup for the region
        for pileup in self.bam.pileup() {
            let pileup = pileup.unwrap();
            let m_pos: u64 = pileup.pos().into();
            // We got to truncate
            if m_pos < start || end < m_pos {
                continue;
            }

            println!("Position: {}, Depth: {}", m_pos, pileup.depth());
            
            // Access additional pileup information as needed
            for alignment in pileup.alignments() {
                // https://docs.rs/rust-htslib/latest/rust_htslib/bam/pileup/struct.Alignment.html
                let m_size = match alignment.indel() {
                    Indel::Ins(size) if size as u64 >= self.params.sizemin => size as i64, // self.make_insertion
                    Indel::Del(size) if size as u64 >= self.params.sizemin => -(size as i64), // self.make_deletion
                    _ => continue,
                };
                println!("indel: {:?}", m_size);
            }
        }
        (Haplotype::blank(3, 0), Haplotype::blank(3, 0))
    }

    // fn hap_deduplicate {}
    // fn read_cluster {}
}

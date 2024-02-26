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

use htslib::{Bam, bam::Read, faidx};
use crate::haplotype::Haplotype;

pub struct BamParser {
    pub bam_name: std::path::PathBuf,
    pub ref_name: std::path::PathBuf,
    bam: mut Bam,
    reference: faidx::Reader,
    params: KDParams,
}

impl BamParser {
    pub fn new (bam_name, ref_name, params) -> Self {
        // Open the BAM or CRAM file
        let mut bam = Bam::from_path(bam_path).unwrap();
        let mut reference = faidx::Reader::from_path(ref_name);
        BamParser {
            bam_name,
            ref_name,
            bam,
            reference,
            params,
        }
    }

    fn find_haps(&self, chrom: String, start: u64, end: u64) -> (Haplotype, Haplotype) {
        // Region of interest (e.g., chromosome:start-end)
        let region = format!("{chrom}:{start}-{end}");

        // Seek to the region of interest
        self.bam.seek(region).unwrap();

        // Iterate over the pileup for the region
        for pileup in bam.pileup() {
            let pileup = pileup.unwrap();
            println!("Position: {}, Depth: {}", pileup.pos(), pileup.depth());
            
            // Access additional pileup information as needed
            for alignment in pileup.alignments() {
                println!("Alignment: {:?}", alignment);
            }
        }

    }

    // fn hap_deduplicate {}
    // fn read_cluster {}

}

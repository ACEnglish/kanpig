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
use crate::kmeans::kmeans;
use crate::kmer::seq_to_kmer;
use crate::pileup::PileupVariant;
use crate::vcf_traits::Svtype;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::faidx;
use std::collections::HashMap;
use std::path::PathBuf;

pub struct BamParser {
    pub bam_name: PathBuf,
    pub ref_name: PathBuf,
    bam: Option<IndexedReader>,
    reference: Option<faidx::Reader>,
    params: KDParams,
    is_open: bool,
}

impl BamParser {
    pub fn new(bam_name: PathBuf, ref_name: PathBuf, params: KDParams) -> Self {
        BamParser {
            bam_name,
            ref_name,
            bam: None,
            reference: None,
            params,
            is_open: false,
        }
    }

    // Open the BAM or CRAM file
    // Put this in a separate method so threads can have their own file handler
    pub fn open(&mut self) {
        self.bam = Some(IndexedReader::from_path(&self.bam_name).unwrap());
        self.reference = Some(faidx::Reader::from_path(&self.ref_name).unwrap());
        self.is_open = true;
    }

    pub fn find_haps(&mut self, chrom: String, start: u64, end: u64) -> (Haplotype, Haplotype) {
        if !self.is_open {
            self.open();
        }
        // We pileup a little outside the region for variants
        let window_start = start - self.params.chunksize;
        let window_end = end + self.params.chunksize;

        let mut tot_cov: u64 = 0;
        if let Err(e) = self.bam.as_mut().expect("Must be opened").fetch((
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
        for pileup in self.bam.as_mut().expect("Must be opened").pileup() {
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

                *p_variants.entry(m_var.clone()).or_insert(0) += 1;
                let qname = alignment.record().qname().to_owned();
                m_reads.entry(qname).or_insert_with(Vec::new).push(m_var);
            }
        }

        let mut m_haps = self.reads_to_haps(m_reads, p_variants, chrom);
        let coverage = tot_cov / (window_end - window_start + (2 * self.params.chunksize));
        if coverage == 0 || m_haps.is_empty() {
            return (
                Haplotype::blank(self.params.kmer, coverage),
                Haplotype::blank(self.params.kmer, coverage),
            );
        };

        let REFTHRESHOLD = 0.85;

        if m_haps.len() == 1 {
            let hap2 = m_haps.pop().unwrap();
            let hap1 = if (hap2.coverage as f64 / coverage as f64) < REFTHRESHOLD {
                Haplotype::blank(self.params.kmer, coverage - hap2.coverage)
            } else {
                hap2.clone()
            };
            return (hap1, hap2);
        }

        self.read_cluster(m_haps, coverage)
    }

    fn reads_to_haps(
        &self,
        reads: HashMap<Vec<u8>, Vec<PileupVariant>>,
        pileups: HashMap<PileupVariant, u64>,
        chrom: String,
    ) -> Vec<Haplotype> {
        let mut hap_parts = HashMap::<PileupVariant, Haplotype>::new();

        for (mut p, _) in pileups.into_iter() {
            // Need to fill in deleted sequence
            if p.indel == Svtype::Del {
                let d_start = p.position;
                let d_end = d_start + p.size.abs() as u64;
                p.sequence = Some(
                    self.reference
                        .as_ref()
                        .unwrap()
                        .fetch_seq(&chrom, d_start as usize, d_end as usize)
                        .unwrap()
                        .to_vec(),
                );
            }

            let n_hap = Haplotype::new(
                seq_to_kmer(&p.sequence.clone().unwrap(), self.params.kmer),
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
            ret.push(cur_hap);
        }

        ret
    }

    /// Cluster multiple haplotypes together to try and reduce them to at most two haplotypes
    fn read_cluster(&self, mut m_haps: Vec<Haplotype>, coverage: u64) -> (Haplotype, Haplotype) {
        let allk = m_haps.iter().map(|i| i.kfeat.clone()).collect::<Vec<_>>();
        let clusts = kmeans(&allk, 2);

        let mut hapA = Vec::<Haplotype>::new();
        let mut hapB = Vec::<Haplotype>::new();
        
        // Sort them so we can pick the highest covered haplotype
        for (idx, hap) in m_haps.drain(..).enumerate() {
            if clusts[0].points_idx.contains(&idx) {
                hapA.push(hap);
            } else {
                hapB.push(hap);
            }
        }
        hapA.sort();
        hapB.sort();

        // Guaranteed to have one, it is going to be the alternate in a 0/1
        let mut hap2 = hapA.pop().unwrap();
        hap2.coverage += hapA.iter().map(|i| i.coverage).sum::<u64>();

        let mut hap1 = if !hapB.is_empty() {
            hapB.pop().unwrap()
            // Now we need to check if its Compound Het or sequencing error should be Hom
        } else {
            // Need to check if its REF & Het or just Hom, if REF, need to set the Coverage correct
            Haplotype::blank(self.params.kmer, 0)
        };
        
        // Consolidate the coverage
        hap1.coverage += hapB.iter().map(|i| i.coverage).sum::<u64>();

        (hap1.clone(), hap2.clone())

        // REF &  HET

        // HOM ALT

        // Compound Het or slightly different Hom?
    }
}

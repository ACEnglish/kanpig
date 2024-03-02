use crate::cli::KDParams;
use crate::haplotype::Haplotype;
use crate::kmeans::kmeans;
use crate::kmer::seq_to_kmer;
use crate::metrics;
use crate::pileup::PileupVariant;
use crate::vcf_traits::Svtype;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::faidx;
use std::collections::HashMap;
use std::path::PathBuf;

pub struct BamParser {
    bam: IndexedReader,
    reference: faidx::Reader,
    params: KDParams,
}

impl BamParser {
    pub fn new(bam_name: PathBuf, ref_name: PathBuf, params: KDParams) -> Self {
        let bam = IndexedReader::from_path(&bam_name).unwrap();
        let reference = faidx::Reader::from_path(&ref_name).unwrap();
        BamParser {
            bam: bam,
            reference: reference,
            params,
        }
    }

    pub fn find_haps(&mut self, chrom: &String, start: u64, end: u64) -> (Haplotype, Haplotype) {
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
                // println!("Pileup {:?}", m_var);
                *p_variants.entry(m_var.clone()).or_insert(0) += 1;
                let qname = alignment.record().qname().to_owned();
                m_reads.entry(qname).or_default().push(m_var);
            }
        }

        // This stuff should be moved to a read_cluster procedure
        // bamparser just creates the pileups and feeds it to something that
        // creates expected haplotypes.
        // The actual cut point is actually somewhere inside reads_to_haps.
        let mut m_haps = self.reads_to_haps(m_reads, p_variants, chrom);
        let coverage = tot_cov / (window_end - window_start);
        // println!("Total coverage: {}", coverage);
        if coverage == 0 || m_haps.is_empty() {
            return (
                Haplotype::blank(self.params.kmer, coverage),
                Haplotype::blank(self.params.kmer, coverage),
            );
        };

        // Nothing to cluster
        if m_haps.len() == 1 {
            let hap2 = m_haps.pop().unwrap();
            let ref_cov = (coverage - hap2.coverage) as f64;
            let ret = match metrics::genotyper(ref_cov, hap2.coverage as f64) {
                // And if Ref, should probably be set to lowq
                metrics::GTstate::Ref | metrics::GTstate::Het => {
                    (Haplotype::blank(self.params.kmer, ref_cov as u64), hap2)
                }
                metrics::GTstate::Hom => (hap2.clone(), hap2),
                _ => panic!("The genotyper can't do this, yet"),
            };
            return ret;
        }

        self.read_cluster(m_haps, coverage)
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
            ret.push(cur_hap);
        }

        // Sort so higher coverage haplotypes are preferred as centroids
        ret.sort_by_key(|i| std::cmp::Reverse(i.coverage));
        ret
    }

    /// Cluster multiple haplotypes together to try and reduce them to at most two haplotypes
    /// This is 'actually' the genotyper. Whatever come out of here is mapped to the variants
    /// So inaccurate descriptions of the two haplotypes can not produce good genotypes.
    fn read_cluster(&self, mut m_haps: Vec<Haplotype>, coverage: u64) -> (Haplotype, Haplotype) {
        let allk = m_haps.iter().map(|i| i.kfeat.clone()).collect::<Vec<_>>();
        let clusts = kmeans(&allk, 2);

        let mut hap_a = Vec::<Haplotype>::new();
        let mut hap_b = Vec::<Haplotype>::new();

        // Sort them so we can pick the highest covered haplotype
        for (idx, hap) in m_haps.drain(..).enumerate() {
            if clusts[0].points_idx.contains(&idx) {
                // println!("Assigning A {:?}", hap);
                hap_a.push(hap);
            } else {
                // println!("Assigning B {:?}", hap);
                hap_b.push(hap);
            }
        }
        hap_a.sort();
        hap_b.sort();

        // Guaranteed to have one cluster
        let mut hap2 = hap_a.pop().unwrap();
        hap2.coverage += hap_a.iter().map(|i| i.coverage).sum::<u64>();

        let mut hap1 = if !hap_b.is_empty() {
            let mut hap_t = hap_b.pop().unwrap();
            hap_t.coverage += hap_b.iter().map(|i| i.coverage).sum::<u64>();
            hap_t
        } else {
            // make a ref haplotype from the remaining coverage
            Haplotype::blank(self.params.kmer, coverage - hap2.coverage)
        };

        // We want the second haplotype to be the higher covered if it is a non-REF
        if hap1.n != 0
            && (hap1.coverage > hap2.coverage
                || (hap1.coverage == hap2.coverage && hap1.n < hap2.n))
        {
            std::mem::swap(&mut hap1, &mut hap2);
        }

        // println!("Working with \n\t{:?}\n\t{:?}", hap1, hap2);
        // println!("Hap coverage: {} - {}:{}", coverage, hap1.coverage, hap2.coverage);

        // First we establish the two possible alt alleles
        // This is a dedup step for when the alt paths are highly similar
        let (hap1, mut hap2) = match metrics::genotyper(hap1.coverage as f64, hap2.coverage as f64)
        {
            // if hap1 == ref, return hap1, hap2. else combine hap2 into hap1 and make return Hap::blank, hap2
            metrics::GTstate::Ref => {
                if hap1.n == 0 {
                    // println!("Think hap2 is a better representative");
                    hap2.coverage += hap1.coverage;
                    (Haplotype::blank(self.params.kmer, 0), hap2)
                } else {
                    // println!("Think hap1 is a better representative");
                    hap1.coverage += hap2.coverage;
                    (Haplotype::blank(self.params.kmer, 0), hap1)
                }
            }
            // combine them (into hap2, the higher covered allele) return hap2 twice
            metrics::GTstate::Hom => {
                if (hap1.size.signum() == hap2.size.signum())
                    && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs())
                        > self.params.sizesim
                {
                    // println!("These should be consolidated");
                    hap2.coverage += hap1.coverage;
                    (Haplotype::blank(self.params.kmer, 0), hap2)
                } else {
                    // println!("Could be a compound het or a fp");
                    (hap1, hap2)
                }
            }
            // If they're highly similar, combine and assume it was a 'noisy' HOM. Otherwise compound het
            metrics::GTstate::Het => {
                if (hap1.size.signum() == hap2.size.signum())
                    && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs())
                        > self.params.sizesim
                {
                    // println!("Turning into Hom Alt");
                    hap2.coverage += hap1.coverage;
                    (Haplotype::blank(self.params.kmer, 0), hap2)
                } else {
                    // println!("Compound Het");
                    (hap1, hap2)
                }
            }
            _ => panic!("The genotyper can't do this, yet"),
        };
        // Now we figure out if the we need two alt alleles or not
        // The reason this takes two steps is the above code is just trying to figure out if
        // there's 1 or 2 alts. Now we figure out if its Het/Hom
        // So essentially we're checking if coverage compared to hap1.coverage/hap2.coverage puts
        // us in ref/het/hom
        let applied_coverage = (hap1.coverage + hap2.coverage) as f64;
        let remaining_coverage = coverage as f64 - applied_coverage;
        match metrics::genotyper(remaining_coverage, applied_coverage) {
            // We need the one higher covered alt
            // and probably should assign this GT as lowq if REF
            metrics::GTstate::Ref | metrics::GTstate::Het => {
                hap2.coverage += hap1.coverage;
                (
                    Haplotype::blank(self.params.kmer, remaining_coverage as u64),
                    hap2,
                )
            }
            metrics::GTstate::Hom => {
                // remaining coverage (if any) is lost here
                if hap1.n == 0 {
                    // HOMALT
                    (hap2.clone(), hap2)
                } else {
                    // Compound Het
                    (hap1, hap2)
                }
            }
            _ => panic!("The genotyper can't do this, yet"),
        }
    }
}

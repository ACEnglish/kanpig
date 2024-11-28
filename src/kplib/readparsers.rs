use crate::kplib::{seq_to_kmer, Haplotype, KDParams, PileupVariant, ReadPileup, Svtype};
use indexmap::{IndexMap, IndexSet};
use rust_htslib::faidx;
use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::{IndexedReader, Read as BamRead},
    tbx::{self, Read as TbxRead},
};
use std::path::PathBuf;

pub type ReadsMap = IndexMap<usize, Vec<usize>>;
pub type PileupSet = IndexSet<PileupVariant>;

pub trait ReadParser {
    fn find_pileups(&mut self, chrom: &str, start: u64, end: u64) -> (Vec<Haplotype>, u64);
}

pub struct BamParser {
    bam: IndexedReader,
    reference: faidx::Reader,
    params: KDParams,
}

impl BamParser {
    pub fn new(
        bam_name: PathBuf,
        ref_name: PathBuf,
        reference: faidx::Reader,
        params: KDParams,
    ) -> Self {
        let mut bam = IndexedReader::from_path(bam_name).unwrap();
        let _ = bam.set_reference(ref_name.clone());
        Self {
            bam,
            reference,
            params,
        }
    }
}

impl ReadParser for BamParser {
    /// Returns all unique haplotypes over a region
    fn find_pileups(&mut self, chrom: &str, start: u64, end: u64) -> (Vec<Haplotype>, u64) {
        // We pileup a little outside the region for variants
        let window_start = start.saturating_sub(self.params.neighdist);
        let window_end = end + self.params.neighdist;

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

        pileups_to_haps(
            chrom,
            reads,
            p_variants,
            coverage,
            &self.reference,
            &self.params,
        )
    }
}

pub struct PlupParser {
    tbx: tbx::Reader,
    reference: faidx::Reader,
    params: KDParams,
}

impl PlupParser {
    /// Creates a new `PlupReader` for a given file path.
    pub fn new(file_path: PathBuf, reference: faidx::Reader, params: KDParams) -> Self {
        let tbx = tbx::Reader::from_path(&file_path).expect("Failed to open TBX file");
        Self {
            tbx,
            reference,
            params,
        }
    }
}

impl ReadParser for PlupParser {
    /// Fetch and parse pileups within a specified genomic interval.
    /// Returns the set of haplotypes
    fn find_pileups(&mut self, chrom: &str, start: u64, end: u64) -> (Vec<Haplotype>, u64) {
        let window_start = start.saturating_sub(self.params.neighdist);
        let window_end = end + self.params.neighdist;

        let tid = self
            .tbx
            .tid(chrom)
            .unwrap_or_else(|_| panic!("Could not resolve '{}' to contig ID", chrom));
        self.tbx
            .fetch(tid, window_start, window_end)
            .expect("Could not fetch region from TBX");

        let mut reads = ReadsMap::new();
        let mut p_variants = PileupSet::new();
        let mut coverage = 0;

        for (qname, line) in self.tbx.records().filter_map(Result::ok).enumerate() {
            if let Some(mut read) = ReadPileup::decode(&line, self.params.sizemin, self.params.sizemax)
            {
                if read.start < window_start && read.end > window_end {
                    coverage += 1;
                    for m_var in read.pileups.drain(..) {
                        if m_var.position >= window_start && m_var.position <= window_end {
                            trace!("{:?}", m_var);
                            let (p_idx, _) = p_variants.insert_full(m_var);
                            reads.entry(qname).or_default().push(p_idx);
                        }
                    }
                }
            }
        }

        pileups_to_haps(
            chrom,
            reads,
            p_variants,
            coverage,
            &self.reference,
            &self.params,
        )
    }
}

pub fn pileups_to_haps(
    chrom: &str,
    reads: ReadsMap,
    mut plups: PileupSet,
    coverage: u64,
    reference: &faidx::Reader,
    params: &KDParams,
) -> (Vec<Haplotype>, u64) {
    let mut hap_parts = Vec::<Haplotype>::with_capacity(plups.len());

    while let Some(mut p) = plups.pop() {
        // Need to fill in deleted sequence
        let sequence = match p.indel {
            Svtype::Del => reference
                .fetch_seq(chrom, p.position as usize, p.end as usize)
                .unwrap()
                .to_vec(),
            Svtype::Ins => p
                .sequence
                .take()
                .expect("Insertions should already have a sequence"),
            _ => panic!("Unknown Svtype"),
        };

        let n_hap = Haplotype::new(
            seq_to_kmer(
                &sequence,
                params.kmer,
                p.indel == Svtype::Del,
                params.maxhom,
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
        let mut cur_hap = Haplotype::blank(params.kmer, 1);
        for p in read_pileups {
            cur_hap.add(&hap_parts[hap_parts.len() - *p - 1]);
        }
        cur_hap.coverage = coverage;
        ret.push(cur_hap);
    }

    ret.sort_by(|a, b| b.cmp(a));
    trace!("{:?}", ret);
    (ret, coverage)
}

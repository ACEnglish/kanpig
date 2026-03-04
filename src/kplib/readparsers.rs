use crate::kplib::{
    pileup::{pileup_finisher, PileupSet, ReadPileup, ReadsMap},
    CoverageTrack, GraphParams, Haplotype, SequenceMeta, VariantGraph,
};
use rust_htslib::{
    bam::ext::BamRecordExtensions,
    bam::{self, IndexedReader, Read as BamRead},
    faidx,
    tbx::{self, Read as TbxRead},
};
use std::path::PathBuf;

pub trait ReadParser {
    /// Pull reads
    fn find_reads(&mut self, chrom: &str, start: u64, end: u64)
        -> (Vec<ReadPileup>, CoverageTrack);
    /// Official Sample Name
    fn get_sample_name(&self) -> String;
    /// Index of the sample - this is for SequenceMeta which holds all samples at once in vectors
    fn get_sample_idx(&self) -> usize;
    /// Needed for initializing the SequenceMeta with the correct array size
    fn get_sample_count(&self) -> usize;
}

pub struct BamParser {
    bam: IndexedReader,
    reference: faidx::Reader,
    sample_name: String,
    sample_idx: usize,
    sample_count: usize,
    params: GraphParams,
}

impl BamParser {
    pub fn new(
        bam_name: PathBuf,
        ref_name: PathBuf,
        reference: faidx::Reader,
        sample_name: String,
        sample_idx: usize,
        sample_count: usize,
        params: GraphParams,
    ) -> Self {
        let mut bam = open_bam(&bam_name).expect("BAM already checked");
        let _ = bam.set_reference(ref_name.clone());
        Self {
            bam,
            reference,
            sample_name,
            sample_idx,
            sample_count,
            params,
        }
    }
}

impl ReadParser for BamParser {
    /// Returns all unique haplotypes over a region
    fn find_reads(
        &mut self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> (Vec<ReadPileup>, CoverageTrack) {
        // We pileup a little outside the region for variants
        let window_start = start.saturating_sub(self.params.neighdist);
        let window_end = end + self.params.neighdist;

        if let Err(e) = self.bam.fetch((&chrom, window_start, window_end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        let seq_meta_template = SequenceMeta::new(self.get_sample_idx(), self.get_sample_count());

        let mut reads: Vec<ReadPileup> = vec![];
        // track the changes made by each read
        // We keep a unique set so we don't have to carry around as much data at first
        // But also so we can reduce the IO for fetching deletions' reference sequence
        // Though, mainly this pattern is a vestage of pre-subinterval work
        let mut read_pileup_lookup = ReadsMap::new();
        let mut p_variants = PileupSet::new();

        let mut coverage: Vec<(u64, u64)> = vec![];

        let mut qname = 0;
        let mut record = bam::Record::new();

        while let Some(r) = self.bam.read(&mut record) {
            r.expect("Failed to parse record");
            if !record.seq().is_empty()
                && record.mapq() >= self.params.mapq
                && (record.flags() & self.params.mapflag) == 0
            {
                coverage.push((
                    record.reference_start() as u64,
                    record.reference_end() as u64,
                ));

                let mut read = ReadPileup::new_record(
                    chrom.to_string(),
                    &record,
                    self.params.sizemin,
                    self.params.sizemax,
                    seq_meta_template.clone(),
                );

                for m_var in read.pileups.drain(..) {
                    if m_var.position >= window_start && m_var.position <= window_end {
                        let (p_idx, _) = p_variants.insert_full(m_var);
                        read_pileup_lookup.entry(qname).or_default().push(p_idx);
                    }
                }

                reads.push(read);
                qname += 1;
            }
        }

        pileup_finisher(
            chrom,
            &mut reads,
            read_pileup_lookup,
            p_variants,
            &self.reference,
        );
        (
            reads,
            CoverageTrack::new(Some(coverage), self.params.neighdist),
        )
    }

    fn get_sample_name(&self) -> String {
        self.sample_name.clone()
    }

    fn get_sample_idx(&self) -> usize {
        self.sample_idx
    }

    fn get_sample_count(&self) -> usize {
        self.sample_count
    }
}

pub struct PlupParser {
    tbx: tbx::Reader,
    reference: faidx::Reader,
    sample_name: String,
    sample_idx: usize,
    sample_count: usize,
    params: GraphParams,
}

impl PlupParser {
    /// Creates a new `PlupReader` for a given file path.
    pub fn new(
        file_path: PathBuf,
        reference: faidx::Reader,
        sample_name: String,
        sample_idx: usize,
        sample_count: usize,
        params: GraphParams,
    ) -> Self {
        let tbx = tbx::Reader::from_path(&file_path).expect("Failed to open TBX file");
        Self {
            tbx,
            reference,
            sample_name,
            sample_idx,
            sample_count,
            params,
        }
    }
}

impl ReadParser for PlupParser {
    /// Fetch and parse pileups within a specified genomic interval.
    /// Returns the set of haplotypes
    fn find_reads(
        &mut self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> (Vec<ReadPileup>, CoverageTrack) {
        let window_start = start.saturating_sub(self.params.neighdist);
        let window_end = end + self.params.neighdist;

        let tid = match self.tbx.tid(chrom) {
            Ok(t) => t,
            Err(_) => return (vec![], CoverageTrack::new(None, 0)),
        };
        self.tbx
            .fetch(tid, window_start, window_end)
            .expect("Could not fetch region from TBX");

        let seq_meta_template = SequenceMeta::new(self.get_sample_idx(), self.get_sample_count());
        let mut reads: Vec<ReadPileup> = vec![];
        let mut read_pileup_lookup = ReadsMap::new();
        let mut p_variants = PileupSet::new();

        let mut coverage: Vec<(u64, u64)> = vec![];

        for (qname, line) in self.tbx.records().filter_map(Result::ok).enumerate() {
            if let Some(mut read) = ReadPileup::new_text(
                &line,
                self.params.sizemin,
                self.params.sizemax,
                seq_meta_template.clone(),
            ) {
                coverage.push((read.start, read.end));

                for m_var in read.pileups.drain(..) {
                    if m_var.position >= window_start && m_var.position <= window_end {
                        let (p_idx, _) = p_variants.insert_full(m_var);
                        read_pileup_lookup.entry(qname).or_default().push(p_idx);
                    }
                }

                reads.push(read);
            }
        }

        pileup_finisher(
            chrom,
            &mut reads,
            read_pileup_lookup,
            p_variants,
            &self.reference,
        );
        (
            reads,
            CoverageTrack::new(Some(coverage), self.params.neighdist),
        )
    }

    fn get_sample_name(&self) -> String {
        self.sample_name.clone()
    }

    fn get_sample_idx(&self) -> usize {
        self.sample_idx
    }

    fn get_sample_count(&self) -> usize {
        self.sample_count
    }
}

/// Open an optionally remote bam file
pub fn open_bam(bam_path: &PathBuf) -> Result<IndexedReader, Box<dyn std::error::Error>> {
    let path_str = bam_path.to_str().ok_or("Invalid UTF-8 in path")?;

    let reader = if path_str.contains("://") {
        let url = url::Url::parse(path_str)?;
        IndexedReader::from_url(&url)?
    } else {
        IndexedReader::from_path(bam_path)?
    };

    Ok(reader)
}

/// Factory function for opening either a bam or a plup
pub fn open_reads(
    reads_path: PathBuf,
    reference_path: PathBuf,
    sample_name: String,
    sample_idx: usize,
    sample_count: usize,
    params: &GraphParams,
) -> Box<dyn ReadParser> {
    let reference = faidx::Reader::from_path(&reference_path).unwrap();
    match reads_path.file_name().and_then(|name| name.to_str()) {
        Some(name) if name.ends_with(".plup.gz") => Box::new(PlupParser::new(
            reads_path,
            reference,
            sample_name,
            sample_idx,
            sample_count,
            params.clone(),
        )),
        _ => Box::new(BamParser::new(
            reads_path,
            reference_path,
            reference,
            sample_name,
            sample_idx,
            sample_count,
            params.clone(),
        )),
    }
}

// Data structure to hold pileup information
// This needs helper methods to subset read information to subintervals so that we can
// create the same information that we need in each of the modes
// So 1. refactor germ mode so it uses these
// Then 2. I won't need to work so hard in trio/mosaic to lift that code over
// I still don't know where that coverge vector is going to come into play
// Its very important for trio/mosaic
#[derive(Clone)]
pub struct ReadData {
    pub reads: Vec<ReadPileup>,
    pub coverages: Vec<CoverageTrack>,
}

impl ReadData {
    ///
    /// Do read trimming and conversion to haplotypes as well as making the
    /// classic coverage vectors. Note that this can only be called on reads
    /// with sequence metadata holding info on a single sample.
    ///
    pub fn subset_to_interval(
        &mut self,
        start: u64,
        end: u64,
        kmer: (u8, u8),
    ) -> (Vec<Haplotype>, Vec<u64>, Vec<usize>) {
        let mut m_haps: Vec<Haplotype> = vec![];
        let mut m_coverage: Vec<u64> = vec![0; self.coverages.len()];
        let mut m_ref_coverage: Vec<usize> = vec![0; self.coverages.len()];

        for read in self.reads.iter() {
            if !read.spans(start, end) {
                continue;
            }
            let sample_idx = read.meta.samples_flag.trailing_zeros() as usize;
            m_coverage[sample_idx] += 1;
            if !read.pileups.is_empty() {
                m_haps.push(Haplotype::from_readpileup(read.trim(start, end), kmer));
            } else {
                m_ref_coverage[sample_idx] += 1
            }
        }

        (m_haps, m_coverage, m_ref_coverage)
    }
}

pub fn collect_read_data(
    samples: &mut Vec<Box<dyn ReadParser>>,
    m_graph: &VariantGraph,
) -> ReadData {
    let mut reads = Vec::<ReadPileup>::new();
    let mut coverages = Vec::<CoverageTrack>::with_capacity(samples.len());
    for samp in samples.iter_mut() {
        let (m_reads, cov_track) = samp.find_reads(&m_graph.chrom, m_graph.start, m_graph.end);
        reads.extend(m_reads);
        coverages.push(cov_track);
    }

    ReadData { reads, coverages }
}

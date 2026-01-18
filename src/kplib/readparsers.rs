use crate::kplib::{
    pileup::{pileup_finisher, ReadPileup, ReadsMap, PileupSet},
    CoverageTrack, GraphParams, Haplotype, HaplotypeMeta, SequenceMeta, VariantGraph
};
use rust_htslib::{
    faidx,
    bam::ext::BamRecordExtensions,
    bam::{self, IndexedReader, Read as BamRead},
    tbx::{self, Read as TbxRead},
};
use std::path::PathBuf;


pub trait ReadParser {
    /// Pull reads
    fn find_reads(&mut self, chrom: &str, start: u64, end: u64) -> (Vec<ReadPileup>, CoverageTrack);
    /// Official Sample Name
    fn get_sample_name(&self) -> String;
    /// Index of the sample - this is for HaplotypeMeta which holds all samples at once in vectors
    fn get_sample_idx(&self) -> usize;
    /// Needed for initializing the HaplotypeMeta with the correct array size
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
    fn find_reads(&mut self, chrom: &str, start: u64, end: u64) -> (Vec<ReadPileup>, CoverageTrack) {
        // We pileup a little outside the region for variants
        let window_start = start.saturating_sub(self.params.neighdist);
        let window_end = end + self.params.neighdist;

        if let Err(e) = self.bam.fetch((&chrom, window_start, window_end)) {
            panic!("Unable to fetch bam {}:{}-{}\n{:?}", chrom, start, end, e)
        };

        let mut seq_meta_template = SequenceMeta::new(self.get_sample_idx(), self.get_sample_count());

        let mut reads : Vec<ReadPileup> = vec![];
        // track the changes made by each read
        // We keep a unique set so we don't have to carry around as much data at first
        // But also so we can reduce the IO for fetching deletions' reference sequence
        // Though, mainly this pattern is a vestage of pre-subinterval work
        let mut read_pileup_lookup = ReadsMap::new();
        let mut p_variants = PileupSet::new();

        let mut coverage : Vec<(u64, u64)> = vec![];

        let mut qname = 0;
        let mut record = bam::Record::new();

        while let Some(r) = self.bam.read(&mut record) {
            r.expect("Failed to parse record");
            if !record.seq().is_empty()
                && record.mapq() >= self.params.mapq
                && (record.flags() & self.params.mapflag) == 0
            {
                coverage.push((record.reference_start() as u64, record.reference_end() as u64));

                let mut read = ReadPileup::new_record(
                    chrom.to_string(),
                    &record,
                    self.params.sizemin,
                    self.params.sizemax,
                    seq_meta_template.clone(),
                );

                // Kind of weird we yoink all the pileups and then place them back later
                // But I guess that's fine
                for m_var in read.pileups.drain(..) {
                    if m_var.position >= window_start && m_var.position <= window_end {
                        let (p_idx, _) = p_variants.insert_full(m_var);
                        read_pileup_lookup.entry(qname).or_default().push(p_idx);
                    }
                }

                qname += 1;
            }
        }
        
        ( 
            pileup_finisher(chrom, reads, read_pileup_lookup, p_variants, &self.reference),
            CoverageTrack::new(Some(coverage), self.params.neighdist)
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
    fn find_reads(&mut self, chrom: &str, start: u64, end: u64) -> (Vec<ReadPileup>, CoverageTrack) {
        let window_start = start.saturating_sub(self.params.neighdist);
        let window_end = end + self.params.neighdist;

        let tid = match self.tbx.tid(chrom) {
            Ok(t) => t,
            Err(_) => return (vec![], CoverageTrack::new(None, 0)),
        };
        self.tbx
            .fetch(tid, window_start, window_end)
            .expect("Could not fetch region from TBX");

        let mut seq_meta_template = SequenceMeta::new(self.get_sample_idx(), self.get_sample_count());
        let mut reads : Vec<ReadPileup> = vec![];
        let mut read_pileup_lookup = ReadsMap::new();
        let mut p_variants = PileupSet::new();

        let mut coverage : Vec<(u64, u64)> = vec![];

        for (qname, line) in self.tbx.records().filter_map(Result::ok).enumerate() {
            if let Some(mut read) =
                ReadPileup::new_text(&line, self.params.sizemin, self.params.sizemax, seq_meta_template.clone())
            {
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

        ( 
            pileup_finisher(chrom, reads, read_pileup_lookup, p_variants, &self.reference),
            CoverageTrack::new(Some(coverage), self.params.neighdist)
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

// Data structure to hold pileup information from multiple samples
#[derive(Clone)]
pub struct ReadData {
    pub reads: Vec<ReadPileup>, // TODO: These should be ReadPileup
    pub coverages: Vec<CoverageTrack>,
}

pub fn collect_read_data(
    samples: &mut Vec<Box<dyn ReadParser>>,
    m_graph: &VariantGraph,
) -> ReadData {
    let mut reads = Vec::<Haplotype>::new();
    let mut coverages = Vec::<CoverageTrack>::with_capacity(samples.len());
    for samp in samples.iter_mut() {
        let (reads, cov_track) = samp.find_reads(&m_graph.chrom, m_graph.start, m_graph.end);
        reads.extend(reads);
        coverages.push(cov_track);
    }

    ReadData {
        reads,
        coverages,
    }
}


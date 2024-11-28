extern crate pretty_env_logger;

use clap::{Parser, Subcommand};

#[derive(Parser, Clone, Debug)]
#[command(name = "kanpig")]
#[command(about = "Kmer ANalysis of PIleups for Genotyping")]
#[command(author = "ACEnglish", version)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

pub trait KanpigParams {
    fn validate(&self) -> bool;
    fn trace(&self) -> bool;
    fn debug(&self) -> bool;
}

#[derive(Subcommand, Debug, Clone)]
pub enum Commands {
    #[command(about = "Genotype SVs")]
    Gt(GTArgs),

    #[command(about = "BAM/CRAM to Pileup Index")]
    Plup(PlupArgs),
}

#[derive(Parser, Debug, Clone)]
pub struct PlupArgs {
    /// Input BAM/CRAM file
    #[arg(short, long)]
    pub bam: std::path::PathBuf,

    /// CRAM file reference
    #[arg(short, long)]
    pub reference: Option<std::path::PathBuf>,

    /// Output file
    #[arg(short, long)]
    pub output: Option<std::path::PathBuf>,

    /// Number of threads
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,

    /// Minimum size of variant to index
    #[arg(long, default_value_t = 50)]
    pub sizemin: u32,

    /// Maximum size of variant to index
    #[arg(long, default_value_t = 10000)]
    pub sizemax: u32,

    /// Minimum mapq of reads to consider
    #[arg(long, default_value_t = 5)]
    pub mapq: u8,

    /// Alignments with flag matching this value are ignored
    #[arg(long, default_value_t = 3840)]
    pub mapflag: u16,

    /// Chunksize in Mbp
    #[arg(long, default_value_t = 25)]
    pub chunk_size: u64,

    /// Verbose logging
    #[arg(long, default_value_t = false)]
    pub debug: bool,

    /// Very Verbose logging
    #[arg(long, default_value_t = false)]
    pub trace: bool,
}

impl KanpigParams for PlupArgs {
    fn trace(&self) -> bool {
        self.trace
    }

    fn debug(&self) -> bool {
        self.debug
    }

    fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.bam.exists() {
            error!("--bam does not exist");
            is_ok = false;
        } else if !self.bam.is_file() {
            error!("--bam is not a file");
            is_ok = false;
        }

        if let Some(path) = &self.reference {
            if !path.exists() {
                error!("--reference does not exist");
                is_ok = false;
            } else if !path.is_file() {
                error!("--reference is not a file");
                is_ok = false;
            }
        }

        if self.sizemin < 20 {
            warn!("--sizemin is recommended to be at least 20");
        }

        is_ok
    }
}

#[derive(Parser, Debug, Clone)]
pub struct GTArgs {
    #[command(flatten)]
    pub io: IOParams,

    #[command(flatten)]
    pub kd: KDParams,
}

#[derive(clap::Args, Clone, Debug)]
pub struct IOParams {
    /// VCF to genotype
    #[arg(short, long)]
    pub input: std::path::PathBuf,

    /// Reads to genotype (.bam, .cram, or .plup.gz)
    #[arg(short, long)]
    pub reads: std::path::PathBuf,

    /// Reference reads are aligned to
    #[arg(short = 'f', long)]
    pub reference: std::path::PathBuf,

    /// Output vcf (unsorted, uncompressed, default stdout)
    #[arg(short, long)]
    pub out: Option<std::path::PathBuf>,

    /// Regions to analyze
    #[arg(long)]
    pub bed: Option<std::path::PathBuf>,

    /// Bed file of non-diploid regions
    #[arg(long)]
    pub ploidy_bed: Option<std::path::PathBuf>,

    /// Sample to apply genotypes to, default first column
    #[arg(long)]
    pub sample: Option<String>,

    /// Number of threads
    #[arg(long, default_value_t = 1)]
    pub threads: usize,

    /// Verbose logging
    #[arg(long, default_value_t = false)]
    pub debug: bool,

    /// Very Verbose logging
    #[arg(long, default_value_t = false)]
    pub trace: bool,
}

#[derive(clap::Args, Clone, Debug)]
pub struct KDParams {
    /// Kmer size for featurization
    #[arg(long, default_value_t = 4)]
    pub kmer: u8,

    /// Minimum distance between variants to create independent graphs
    #[arg(long, default_value_t = 1000)]
    pub neighdist: u64,

    /// Only analyze variants with PASS FILTER
    #[arg(long, default_value_t = false)]
    pub passonly: bool,

    /// Minimum size of variant to analyze
    #[arg(long, default_value_t = 50)]
    pub sizemin: u64,

    /// Maximum size of variant to analyze
    #[arg(long, default_value_t = 10000)]
    pub sizemax: u64,

    /// Maximum number of paths in a graph to traverse
    #[arg(long, default_value_t = 5000)]
    pub maxpaths: u64,

    /// Minimum sequence similarity for paths
    #[arg(long, default_value_t = 0.90)]
    pub seqsim: f32,

    /// Minimum size similarity for paths
    #[arg(long, default_value_t = 0.90)]
    pub sizesim: f32,

    /// Minimum frequency of kmer
    #[arg(long, default_value_t = 2)]
    pub minkfreq: u64,

    /// Haplotype size similarity collapse threshold (off=1)
    #[arg(long, default_value_t = 1.0)]
    pub hapsim: f32,

    /// Scoring penalty for 'gaps'
    #[arg(long, default_value_t = 0.02)]
    pub gpenalty: f32,

    /// Scoring penalty for 'fns'
    #[arg(long, default_value_t = 0.10)]
    pub fpenalty: f32,

    /// Maximum number of FNs allowed in a chunk
    #[arg(long, default_value_t = 3)]
    pub fnmax: usize,

    /// Maximum number of pileups in a chunk to attempt partials
    #[arg(long, default_value_t = 100)]
    pub pileupmax: usize,

    /// Search for a 1-to-1 match before graph traversal
    #[arg(long, default_value_t = false)]
    pub try_exact: bool,

    /// Prune paths which don't traverse 1-to-1 nodes
    #[arg(long, default_value_t = false)]
    pub prune: bool,

    /// Minimum mapq of reads to consider
    #[arg(long, default_value_t = 5)]
    pub mapq: u8,

    /// Alignments with flag matching this value are ignored
    #[arg(long, default_value_t = 3840)]
    pub mapflag: u16,

    /// Maximum homopolymer length to kmerize (off=0)
    #[arg(long, default_value_t = 0)]
    pub maxhom: usize,
}

impl KanpigParams for GTArgs {
    fn trace(&self) -> bool {
        self.io.trace
    }

    fn debug(&self) -> bool {
        self.io.debug
    }
    /// Validate command line arguments
    fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.io.input.exists() {
            error!("--input does not exist");
            is_ok = false;
        } else if !self.io.input.is_file() {
            error!("--input is not a file");
            is_ok = false;
        }

        if !self.io.reads.exists() {
            error!("--reads does not exist");
            is_ok = false;
        } else if !self.io.reads.is_file() {
            error!("--reads is not a file");
            is_ok = false;
        }

        if !self.io.reference.exists() {
            error!("--reference does not exist");
            is_ok = false;
        } else if !self.io.reference.is_file() {
            error!("--reference is not a file");
            is_ok = false;
        }

        let mut fai_path = self.io.reference.clone();
        fai_path.set_file_name(format!(
            "{}.fai",
            fai_path.file_name().unwrap().to_string_lossy()
        ));
        if !fai_path.exists() {
            error!("--reference index (.fai) does not exist");
            is_ok = false;
        } else if !fai_path.is_file() {
            error!("--reference index is not a file");
            is_ok = false;
        }

        if let Some(bed_file) = &self.io.bed {
            if !bed_file.exists() {
                error!("--bed does not exist");
                is_ok = false;
            } else if !bed_file.is_file() {
                error!("--bed is not a file");
                is_ok = false;
            }
        }

        if self.kd.sizemin < 20 {
            warn!("--sizemin is recommended to be at least 20");
        }

        if self.kd.kmer >= 8 {
            warn!("--kmer above 8 becomes memory intensive");
        }

        if self.kd.kmer < 1 {
            error!("--kmer must be at least 1");
            is_ok = false;
        }

        if self.kd.sizesim < 0.0 || self.kd.sizesim > 1.0 {
            error!("--sizesim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.kd.seqsim < 0.0 || self.kd.seqsim > 1.0 {
            error!("--seqsim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.kd.hapsim < 0.0 || self.kd.hapsim > 1.0 {
            error!("--hapsim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.kd.maxpaths < 1 {
            error!("--maxpaths must be at least 1");
            is_ok = false;
        }

        if self.io.threads < 1 {
            error!("--threads must be at least 1");
            is_ok = false;
        }

        is_ok
    }
}

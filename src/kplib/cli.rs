extern crate pretty_env_logger;

use clap::Parser;

#[derive(Parser, Clone, Debug)]
#[command(author = "ACEnglish", version)]
pub struct ArgParser {
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

    /// Reads to genotype
    #[arg(short, long)]
    pub bam: std::path::PathBuf,

    /// Reference bam is aligned to
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
}

#[derive(clap::Args, Clone, Debug)]
pub struct KDParams {
    /// Kmer size for featurization
    #[arg(long, default_value_t = 4)]
    pub kmer: u8,

    /// Minimum distance between variants to create independent graphs
    #[arg(long, default_value_t = 1000)]
    pub chunksize: u64,

    /// Only analyze reads with PASS FILTER
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
    #[arg(long, default_value_t = 1)]
    pub minkfreq: u64,

    /// Haplotype size similarity collapse threshold (off=1)
    #[arg(long, default_value_t = 1.0)]
    pub hapsim: f32,

    /// Scoring penalty for 'gaps'
    #[arg(long, default_value_t = 0.01)]
    pub gpenalty: f32,

    /// Scoring penalty for 'fns'
    #[arg(long, default_value_t = 0.10)]
    pub fpenalty: f32,

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

    /// Don't require alignments to span vargraph region
    #[arg(long, default_value_t = true)]
    pub spanoff: bool,

    /// Maximum homopolymer length to kmerize (off=0)
    #[arg(long, default_value_t = 0)]
    pub maxhom: usize,
}

impl ArgParser {
    /// Validate command line arguments
    pub fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.io.input.exists() {
            error!("--input does not exist");
            is_ok = false;
        } else if !self.io.input.is_file() {
            error!("--input is not a file");
            is_ok = false;
        }

        if !self.io.bam.exists() {
            error!("--bam does not exist");
            is_ok = false;
        } else if !self.io.bam.is_file() {
            error!("--bam is not a file");
            is_ok = false;
        }

        if !self.io.reference.exists() {
            error!("--reference does not exist");
            is_ok = false;
        } else if !self.io.reference.is_file() {
            error!("--reference is not a file");
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

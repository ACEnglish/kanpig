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

    /// Output vcf (unsorted)
    #[arg(short, long)]
    pub out: std::path::PathBuf,

    /// Regions to analyze
    #[arg(long)]
    pub bed: Option<std::path::PathBuf>,

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
    #[arg(long, default_value_t = 100)]
    pub chunksize: u64,

    /// Only analyze reads with PASS FILTER
    #[arg(long, default_value_t = false)]
    pub passonly: bool,

    /// Minimum size of variant to analyze
    #[arg(long, default_value_t = 50)]
    pub sizemin: u64,

    // change this to 10k? 
    /// Maximum size of variant to analyze
    #[arg(long, default_value_t = 50000)]
    pub sizemax: u64,

    /// Maximum number of paths in a graph to traverse
    #[arg(long, default_value_t = 1000)]
    pub maxpaths: u64,

    /// Minimum sequence similarity for paths
    #[arg(long, default_value_t = 0.90)]
    pub seqsim: f32,

    /// Minimum size similarity for paths
    #[arg(long, default_value_t = 0.95)]
    pub sizesim: f32,

    /// Minimum frequency of kmer
    #[arg(long, default_value_t = 1)]
    pub minkfreq: u64,
    
    /// Haplotype size similarity collapse threshold
    #[arg(long, default_value_t = 0.95)]
    pub hapsim: f32,

    /// Search for a 1-to-1 match before graph traversal
    #[arg(long, default_value_t = false)]
    pub try_exact: bool,

    /// Don't prune paths which don't traverse 1-to-1 nodes
    #[arg(long, default_value_t = false)]
    pub no_prune: bool,

    /// Minimum mapq of reads to consider
    #[arg(long, default_value_t = 5)]
    pub mapq: u8,

    /// Alignments with flag matching this value are ignored
    #[arg(long, default_value_t = 3840)]
    pub mapflag: u16,
}

impl ArgParser {
    /// Validate command line arguments
    pub fn validate(&self) -> bool {
        let mut is_ok = true;

        if !self.io.input.exists() {
            error!("--input does not exist");
            is_ok = false;
        }

        if !self.io.bam.exists() {
            error!("--bam does not exist");
            is_ok = false;
        }

        if !self.io.reference.exists() {
            error!("--reference does not exist");
            is_ok = false;
        }

        if let Some(bed_file) = &self.io.bed {
            if !bed_file.exists() {
                error!("--bed does not exist");
                is_ok = false;
            }
        }

        if self.kd.sizemin < 20 {
            warn!("--sizemin is recommended to be at least 20");
        }

        if self.kd.kmer >= 8 {
            warn!("--kmer above 8 becomes memory intensive");
        }

        is_ok
    }
}

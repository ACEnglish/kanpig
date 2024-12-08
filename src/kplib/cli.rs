extern crate pretty_env_logger;

use clap::{Parser, Subcommand};
use rust_htslib::tbx::{self, Read as TbxRead};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Parser, Clone, Debug)]
#[command(name = "kanpig")]
#[command(about = "Kmer ANalysis of PIleups for Genotyping")]
#[command(author = "ACEnglish", version)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

pub trait KanpigParams: std::fmt::Debug {
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

#[derive(Parser, Serialize, Deserialize, Debug, Clone)]
pub struct PlupArgs {
    /// Input BAM/CRAM file
    #[arg(short, long)]
    pub bam: PathBuf,

    /// Reference file for CRAMs
    #[arg(short, long)]
    pub reference: Option<PathBuf>,

    /// Output plup (unsorted, uncompressed, default stdout)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

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

        is_ok &= validate_file(&self.bam, "--bam");

        if let Some(ref_path) = &self.reference {
            is_ok &= validate_reference(ref_path);
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
    pub input: PathBuf,

    /// Reads to genotype (indexed .bam, .cram, or .plup.gz)
    #[arg(short, long)]
    pub reads: PathBuf,

    /// Reference reads are aligned to
    #[arg(short = 'f', long)]
    pub reference: PathBuf,

    /// Output vcf (unsorted, uncompressed, default stdout)
    #[arg(short, long)]
    pub out: Option<PathBuf>,

    /// Number of threads
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,

    /// Regions to analyze
    #[arg(long)]
    pub bed: Option<PathBuf>,

    /// Bed file of non-diploid regions
    #[arg(long)]
    pub ploidy_bed: Option<PathBuf>,

    /// Sample to apply genotypes to, default first column
    #[arg(long)]
    pub sample: Option<String>,

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
    pub sizemin: u32,

    /// Maximum size of variant to analyze
    #[arg(long, default_value_t = 10000)]
    pub sizemax: u32,

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

    /// Read clustering weight for haplotagged reads (off=0, full=1)
    #[arg(long, default_value_t = 1.0)]
    pub hps_weight: f32,

    /// Minimum mapq of reads to consider
    #[arg(long, default_value_t = 5)]
    pub mapq: u8,

    /// Alignments with flag matching this value are ignored
    #[arg(long, default_value_t = 3840)]
    pub mapflag: u16,

    /// (Experimental) Search for a 1-to-1 match before graph traversal
    #[arg(long, default_value_t = false)]
    pub try_exact: bool,

    /// (Experimental) Prune paths which don't traverse 1-to-1 nodes
    #[arg(long, default_value_t = false)]
    pub prune: bool,

    /// (Experimental) Maximum homopolymer length to kmerize (off=0)
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

        is_ok &= validate_file(&self.io.input, "--input");
        is_ok &= validate_reads(&self.io.reads, self);
        is_ok &= validate_reference(&self.io.reference);

        if let Some(bed_file) = &self.io.bed {
            is_ok &= validate_file(bed_file, "--bed");
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

/// Helper function to validate a file's existence and type
fn validate_file(path: &Path, label: &str) -> bool {
    if !path.exists() {
        error!("{} does not exist", label);
        return false;
    }
    if !path.is_file() {
        error!("{} is not a file", label);
        return false;
    }
    true
}

/// Helper function to validate reads (.bam, .cram, or .plup.gz)
fn validate_reads(reads: &Path, params: &GTArgs) -> bool {
    let mut is_ok = validate_file(reads, "--reads");

    let file_path = reads.to_str().unwrap_or_default();
    if file_path.ends_with(".bam") || file_path.ends_with(".cram") {
        let index_extensions = [".bai", ".crai"];
        let index_exists = index_extensions.iter().any(|ext| {
            let index_path = format!("{}{}", file_path, ext);
            Path::new(&index_path).exists()
        });

        if !index_exists {
            error!(
                "Index file for {} does not exist. Expected one of: {}",
                file_path,
                index_extensions.join(", ")
            );
            is_ok = false;
        }
    } else if file_path.ends_with(".plup.gz") {
        let tbi_path = format!("{}.tbi", file_path);
        if !Path::new(&tbi_path).exists() {
            error!("Index file {} does not exist", tbi_path);
            is_ok = false;
        } else {
            let tbx = tbx::Reader::from_path(file_path).expect("Failed to open TBX file");
            let header = tbx.header();
            if header.len() != 1 {
                error!("Malformed plup.gz header. Unable to validate parameters");
            } else {
                match serde_json::from_str::<PlupArgs>(&header[0][2..]) {
                    Ok(plup_args) => {
                        if plup_args.sizemin != params.kd.sizemin {
                            warn!(
                                "--reads created with plup --sizemin {} != gt --sizemin {}",
                                plup_args.sizemin, params.kd.sizemin
                            );
                        }

                        if plup_args.sizemax != params.kd.sizemax {
                            warn!(
                                "--reads created with plup --sizemax {} != gt --sizemax {}",
                                plup_args.sizemax, params.kd.sizemax
                            );
                        }

                        if plup_args.mapq != params.kd.mapq {
                            warn!(
                                "--reads created with plup --mapq {} != gt --mapq {}",
                                plup_args.mapq, params.kd.mapq
                            );
                        }

                        if plup_args.mapflag != params.kd.mapflag {
                            warn!(
                                "--reads created plup --mapflag {} != gt --mapflag {}",
                                plup_args.mapflag, params.kd.mapflag
                            );
                        }
                    }
                    Err(e) => {
                        error!(
                            "Failed to parse plup.gz header for parameter validation: {}",
                            e
                        );
                    }
                }
            }
        }
    } else {
        error!("Unsupported file type: {}", file_path);
        is_ok = false;
    }

    is_ok
}

/// Checks reference and its .fai index
fn validate_reference(reference: &Path) -> bool {
    let mut is_ok = validate_file(reference, "--reference");

    let mut fai_path = reference.to_path_buf();
    fai_path.set_file_name(format!(
        "{}.fai",
        fai_path.file_name().unwrap().to_string_lossy()
    ));
    is_ok &= validate_file(&fai_path, "--reference index (.fai)");

    is_ok
}

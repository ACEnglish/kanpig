extern crate pretty_env_logger;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(author = "ACEnglish", version)]
pub struct ArgParser {
    #[command(flatten)]
    pub io: IOParams,

    #[command(flatten)]
    pub kd: KDParams,
}

#[derive(clap::Args, Debug)]
pub struct IOParams {
    #[arg(short, long)]
    pub input: std::path::PathBuf,

    #[arg(long)]
    pub vcf: Option<std::path::PathBuf>,

    #[arg(long)]
    pub bam: Option<std::path::PathBuf>,

    #[arg(short = 'f', long)]
    pub reference: Option<std::path::PathBuf>,

    #[arg(short, long)]
    pub out: std::path::PathBuf,

    #[arg(short, long)]
    pub regions: Option<std::path::PathBuf>,

    #[arg(short, long)]
    pub sample: Option<String>,
}

#[derive(clap::Args, Debug, Clone)]
pub struct KDParams {
    #[arg(long, default_value_t = 4)]
    pub kmer: u8,

    #[arg(long, default_value_t = false)]
    pub passonly: bool,

    #[arg(long, default_value_t = 20)]
    pub sizemin: u64,

    #[arg(long, default_value_t = 50000)]
    pub sizemax: u64,

    #[arg(long, default_value_t = 1000)]
    pub maxpaths: usize,

    #[arg(long, default_value_t = 0.90)]
    pub cossim: f32,

    #[arg(long, default_value_t = 0.90)]
    pub pctsize: f32,

    #[arg(long, default_value_t = 2000)]
    pub wcoslen: usize,

    #[arg(long, default_value_t = 5)]
    pub n_tries: usize,

    #[arg(long, default_value_t = 100)]
    pub chunksize: usize,
}

impl ArgParser {
    /// Validate command line arguments
    pub fn validate(&self) -> bool {
        let mut is_ok = true;

        if self.io.bam.is_some() == self.io.vcf.is_some() {
            error!("One of --vcf xor --bam must be provided");
            is_ok = false;
        }

        if self.io.bam.is_some() && self.io.reference.is_none() {
            error!("--bam requires --reference");
            is_ok = false;
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

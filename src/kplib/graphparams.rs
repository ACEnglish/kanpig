///
/// Parameters relevant for variants, reads, and graph that are used by (most) commands
///
#[derive(clap::Args, Clone, Debug)]
pub struct GraphParams {
    /// Only analyze variants with PASS FILTER
    #[arg(long, default_value_t = false, help_heading = "Variants & Reads")]
    pub passonly: bool,

    /// Maximum variant distance within graphs
    #[arg(long, default_value_t = 1000, help_heading = "Variants & Reads")]
    pub neighdist: u64,

    /// Minimum size of variant to analyze
    #[arg(long, default_value_t = 50, help_heading = "Variants & Reads")]
    pub sizemin: u32,

    /// Maximum size of variant to analyze
    #[arg(long, default_value_t = 10000, help_heading = "Variants & Reads")]
    pub sizemax: u32,

    /// Minimum mapq score for reads
    #[arg(long, default_value_t = 5, help_heading = "Variants & Reads")]
    pub mapq: u8,

    /// Ignore alignments matching flag
    #[arg(long, default_value_t = 3840, help_heading = "Variants & Reads")]
    pub mapflag: u16,

    /// Minimum sequence similarity for paths
    #[arg(long, default_value_t = 0.90, help_heading = "Scoring / Advanced")]
    pub seqsim: f32,

    /// Minimum size similarity for paths
    #[arg(long, default_value_t = 0.90, help_heading = "Scoring / Advanced")]
    pub sizesim: f32,

    /// Scoring penalty for gaps
    #[arg(long, default_value_t = 0.02, help_heading = "Scoring / Advanced")]
    pub gpenalty: f32,

    /// Scoring penalty for FNs
    #[arg(long, default_value_t = 0.10, help_heading = "Scoring / Advanced")]
    pub fpenalty: f32,

    /// Kmer size for featurization
    #[arg(long, default_value_t = 4, help_heading = "Scoring / Advanced")]
    pub kmer: u8,

    /// Minimum frequency of kmers
    #[arg(long, default_value_t = 2, help_heading = "Scoring / Advanced")]
    pub minkfreq: u64,

    /// Maximum graph size to search; otherwise perform 1-to-1
    #[arg(long, default_value_t = 5000, help_heading = "Scoring / Advanced")]
    pub maxnodes: usize,

    /// Maximum paths to traverse per graph
    #[arg(long, default_value_t = 5000, help_heading = "Scoring / Advanced")]
    pub maxpaths: u64,

    /// Maximum pileups allowed for partials matching
    #[arg(long, default_value_t = 100, help_heading = "Scoring / Advanced")]
    pub pileupmax: usize,

    /// Maximum FNs allowed in a path
    #[arg(long, default_value_t = 3, help_heading = "Scoring / Advanced")]
    pub fnmax: usize,

    /// Prefer simplier paths during scoring
    #[arg(long, default_value_t = false, help_heading = "Scoring / Advanced")]
    pub squish: bool,

    /// (Experimental) Restrict to 1-to-1 haplotype/node matching
    #[arg(long, default_value_t = false, help_heading = "Scoring / Advanced")]
    pub one_to_one: bool,

    /// (Experimental) Limit homopolymer length (off=0)
    #[arg(long, default_value_t = 0, help_heading = "Scoring / Advanced")]
    pub maxhom: usize,
}

impl Default for GraphParams {
    fn default() -> Self {
        Self {
            passonly: false,
            neighdist: 1000,
            sizemin: 50,
            sizemax: 10000,
            mapq: 5,
            mapflag: 3840,
            seqsim: 0.90,
            sizesim: 0.90,
            gpenalty: 0.02,
            fpenalty: 0.10,
            kmer: 4,
            minkfreq: 2,
            maxnodes: 5000,
            maxpaths: 5000,
            pileupmax: 100,
            fnmax: 3,
            squish: false,
            one_to_one: false,
            maxhom: 0,
        }
    }
}

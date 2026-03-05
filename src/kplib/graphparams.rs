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

    /// Minimum kmer similarity for paths
    #[arg(long, default_value_t = 0.60, help_heading = "Graph")]
    pub seqsim: f32,

    /// Minimum size similarity for paths
    #[arg(long, default_value_t = 0.90, help_heading = "Graph")]
    pub sizesim: f32,

    /// Scoring penalty for gaps
    #[arg(long, default_value_t = 0.02, help_heading = "Graph")]
    pub gpenalty: f32,

    /// Scoring penalty for FNs
    #[arg(long, default_value_t = 0.10, help_heading = "Graph")]
    pub fpenalty: f32,

    /// Kmer size for coarse comparisons (max 32)
    #[arg(long, default_value_t = 16, help_heading = "Graph")]
    pub coarse_kmer: u8,

    /// Kmer size for fine comparisons (max 32)
    #[arg(long, default_value_t = 4, help_heading = "Graph")]
    pub fine_kmer: u8,

    /// Maximum graph size to search; otherwise perform 1-to-1
    #[arg(long, default_value_t = 5000, help_heading = "Graph")]
    pub maxnodes: usize,

    /// Maximum paths to traverse per graph
    #[arg(long, default_value_t = 5000, help_heading = "Graph")]
    pub maxpaths: u64,

    /// Maximum pileups allowed for partials matching
    #[arg(long, default_value_t = 100, help_heading = "Graph")]
    pub pileupmax: usize,

    /// Maximum FNs allowed in a path
    #[arg(long, default_value_t = 3, help_heading = "Graph")]
    pub fnmax: usize,

    /// Prefer simplier paths during scoring
    #[arg(long, default_value_t = false, help_heading = "Graph")]
    pub squish: bool,

    /// Restrict to 1-to-1 haplotype/node matching
    #[arg(long, default_value_t = false, help_heading = "Graph")]
    pub one_to_one: bool,

    /// Minimum coverage to attempt building haplotypes
    #[arg(long, default_value_t = 1, help_heading = "Graph")]
    pub mincoverage: usize,

    /// Maximum coverage to attempt building haplotypes
    #[arg(long, default_value_t = 1000, help_heading = "Graph")]
    pub maxcoverage: usize,

    /// [Experimental] Subset graphs to subintervals around pileups
    #[arg(long, default_value_t = false, help_heading = "Graph")]
    pub subintv: bool,
}

impl GraphParams {
    pub fn validate(&self) -> bool {
        let mut is_ok = true;
        if self.sizemin < 10 {
            warn!("--sizemin is recommended to be at least 10");
        }

        if self.coarse_kmer > 32 || self.fine_kmer > 32 {
            error!("--kmer must be below 32");
            is_ok = false;
        }

        if self.coarse_kmer < 1 || self.fine_kmer < 1 {
            error!("--kmer must be at least 1");
            is_ok = false;
        }

        if self.sizemin < self.coarse_kmer.min(self.fine_kmer).into() {
            error!("--sizemin must be ≥ --kmer");
            is_ok = false;
        }

        if self.sizesim < 0.0 || self.sizesim > 1.0 {
            error!("--sizesim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.seqsim < 0.0 || self.seqsim > 1.0 {
            error!("--seqsim must be between 0.0 and 1.0");
            is_ok = false;
        }

        if self.maxpaths < 1 {
            error!("--maxpaths must be at least 1");
            is_ok = false;
        }

        is_ok
    }
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
            seqsim: 0.80,
            sizesim: 0.85,
            gpenalty: 0.02,
            fpenalty: 0.10,
            coarse_kmer: 16,
            fine_kmer: 4,
            maxnodes: 5000,
            maxpaths: 5000,
            pileupmax: 100,
            fnmax: 3,
            squish: false,
            one_to_one: false,
            mincoverage: 1,
            maxcoverage: 1000,
            subintv: false,
        }
    }
}

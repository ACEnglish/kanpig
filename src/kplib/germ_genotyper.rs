use ordered_float::OrderedFloat;
use serde::Deserialize;
use statrs::{
    distribution::{Binomial, Discrete},
    function::gamma::ln_gamma,
};
use std::fs;
use std::path::PathBuf;

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub enum GTstate {
    Ref,
    Het,
    Hom,
    Non,
    //Hemi should be a thing
}

pub struct GenotypeResult {
    pub state: GTstate,
    pub gq: f64,
    pub sq: f64,
}

#[derive(Debug, Clone, Copy, Deserialize)]
pub enum GenotypeMode {
    Beta,
    Bino,
    Phased,
}

impl GenotypeMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            GenotypeMode::Beta => "Beta",
            GenotypeMode::Bino => "Bino",
            GenotypeMode::Phased => "Phased",
        }
    }

    pub fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "Beta" => Ok(GenotypeMode::Beta),
            "Bino" => Ok(GenotypeMode::Bino),
            "Phased" => Ok(GenotypeMode::Phased),
            other => Err(format!(
                "Invalid GenotypeMode '{}'; expected 'Beta', 'Bino', or 'Phased'",
                other
            )),
        }
    }
}

#[derive(Debug, Deserialize, Clone)]
pub struct GenotyperConfig {
    pub mode: GenotypeMode,
    pub mixture_fractions: Vec<f64>,
    pub means: Vec<f64>,
    pub precisions: Vec<f64>,
    pub calibration_table: Vec<(f64, f64)>,
}

impl GenotyperConfig {
    /// Creates a new Genotyper from a JSON config file
    ///
    /// # Parameters
    /// - `config_path`: Path to the JSON configuration file
    ///
    /// # Returns
    /// Result containing the GenotyperConfig or an error if the file cannot be read/parsed
    pub fn from_config_file(config_path: PathBuf) -> Result<Self, Box<dyn std::error::Error>> {
        let config_str = fs::read_to_string(config_path)?;
        let config: GenotyperConfig = serde_json::from_str(&config_str)?;

        if config.mixture_fractions.len() != 3
            || config.means.len() != 3
            || config.precisions.len() != 3
        {
            return Err("Config must contain exactly 3 values for each parameter".into());
        }

        Ok(config)
    }

    /// Creates a new Genotyper, optionally loading from a config file
    ///
    /// # Parameters
    /// - `config_path`: Optional path to a JSON configuration file
    ///
    /// # Returns
    /// A Genotyper instance (falls back to defaults if path is None or loading fails)
    pub fn from_optional_config(config_path: Option<PathBuf>) -> Self {
        match config_path {
            Some(path) => Self::from_config_file(path).unwrap_or_else(|e| {
                error!("Warning: Failed to load config ({}), using defaults", e);
                Self::default()
            }),
            None => Self::default(),
        }
    }
}

impl Default for GenotyperConfig {
    fn default() -> Self {
        Self {
            mode: GenotypeMode::Beta,
            mixture_fractions: vec![0.33, 0.34, 0.33],
            means: vec![0.03, 0.50, 0.97],
            precisions: vec![25.0, 5.0, 10.0],
            calibration_table: vec![],
        }
    }
}

pub struct Genotyper {
    pub config: GenotyperConfig,
}

impl Genotyper {
    // From config
    pub fn from_config_file(config_path: Option<PathBuf>) -> Self {
        Self {
            config: GenotyperConfig::from_optional_config(config_path),
        }
    }

    pub fn from_config(config: GenotyperConfig) -> Self {
        Self { config }
    }

    /// Determines the genotype state based on coverage values for two alternate alleles.
    /// The genotype state can be one of three: reference (Ref), heterozygous (Het), or homozygous (Hom).
    /// If both coverage values are zero, the state is `Non`.
    ///
    /// # Parameters
    /// - `ref_cov`: The coverage value for the reference allele.
    /// - `alt_cov`: The coverage value for the alternate allele.
    ///
    /// # Returns
    /// A `GenotypeResult` containing the genotype state, GQ, and SQ
    ///
    /// # Panics
    /// This function will panic if an invalid state is encountered, which should be impossible under normal circumstances.
    pub fn genotype(&self, ref_cov: u64, alt_cov1: u64, alt_cov2: u64) -> GenotypeResult {
        let alt_cov = alt_cov1 + alt_cov2;
        let tot_cov = ref_cov + alt_cov1 + alt_cov2;
        if tot_cov == 0 {
            return GenotypeResult {
                state: GTstate::Non,
                gq: 0.0,
                sq: 0.0,
            };
        }

        let scores = match self.config.mode {
            GenotypeMode::Beta => self.beta_genotype_scores(ref_cov, alt_cov),
            GenotypeMode::Bino => self.bino_genotype_scores(ref_cov, alt_cov),
            GenotypeMode::Phased => self.phased_genotype_scores(ref_cov, alt_cov1, alt_cov2),
        };

        let state = match scores
            .iter()
            .enumerate()
            .max_by_key(|&(_, &x)| OrderedFloat(x))
            .map(|(i, _)| i)
        {
            Some(0) => GTstate::Ref,
            Some(1) => GTstate::Het,
            Some(2) => GTstate::Hom,
            _ => panic!("not possible"),
        };
        let (gq, sq) = self.genotype_quals(scores);
        GenotypeResult { state, gq, sq }
    }

    fn beta_binomial_ln_pmf(&self, k: u64, n: u64, alpha: f64, beta: f64) -> f64 {
        let k = k as f64;
        let n = n as f64;

        // log C(n, k) = log(n!) - log(k!) - log((n-k)!)
        let log_binom_coef = ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0);

        // log B(k+α, n-k+β) = log Γ(k+α) + log Γ(n-k+β) - log Γ(n+α+β)
        let log_beta_num =
            ln_gamma(k + alpha) + ln_gamma(n - k + beta) - ln_gamma(n + alpha + beta);

        // log B(α, β) = log Γ(α) + log Γ(β) - log Γ(α+β)
        let log_beta_denom = ln_gamma(alpha) + ln_gamma(beta) - ln_gamma(alpha + beta);

        log_binom_coef + log_beta_num - log_beta_denom
    }

    fn bino_genotype_scores(&self, ref_cov: u64, alt_cov: u64) -> [f64; 3] {
        let error_rate = 0.03;

        // Prior probabilities
        let prior_homref = self.config.mixture_fractions[0].ln();
        let prior_het = self.config.mixture_fractions[1].ln();
        let prior_homalt = self.config.mixture_fractions[2].ln();

        let n = ref_cov + alt_cov;
        if n == 0 {
            return [0.0, 0.0, 0.0];
        }

        // Create binomial distributions for each genotype
        let binom_homref = Binomial::new(error_rate, n).unwrap();
        let binom_het = Binomial::new(0.5, n).unwrap();
        let binom_homalt = Binomial::new(0.98, n).unwrap();

        // Calculate log-likelihoods
        let ll_homref = binom_homref.ln_pmf(alt_cov);
        let ll_het = binom_het.ln_pmf(alt_cov);
        let ll_homalt = binom_homalt.ln_pmf(alt_cov);

        // Posterior = Prior + Likelihood
        [
            prior_homref + ll_homref,
            prior_het + ll_het,
            prior_homalt + ll_homalt,
        ]
    }

    fn beta_genotype_scores(&self, ref_cov: u64, alt_cov: u64) -> [f64; 3] {
        let total = ref_cov + alt_cov;

        if total == 0 {
            return [0.0, 0.0, 0.0]; // keep flat prior for missing GT to keep previous GQ < 5 threshold for LOWGQ filter
        }

        let alpha: Vec<f64> = self
            .config
            .means
            .iter()
            .zip(self.config.precisions.iter())
            .map(|(m, n)| m * n)
            .collect();
        let beta: Vec<f64> = self
            .config
            .means
            .iter()
            .zip(self.config.precisions.iter())
            .map(|(m, n)| (1.0 - m) * n)
            .collect();

        [
            self.config.mixture_fractions[0].ln()
                + self.beta_binomial_ln_pmf(alt_cov, total, alpha[0], beta[0]),
            self.config.mixture_fractions[1].ln()
                + self.beta_binomial_ln_pmf(alt_cov, total, alpha[1], beta[1]),
            self.config.mixture_fractions[2].ln()
                + self.beta_binomial_ln_pmf(alt_cov, total, alpha[2], beta[2]),
        ]
    }

    fn phased_genotype_scores(
        &self,
        ref_reads: u64,     // Reads supporting reference allele
        allele1_reads: u64, // Reads supporting allele1 (haplotype 1)
        allele2_reads: u64, // Reads supporting allele2 (haplotype 2)
    ) -> [f64; 3] {
        let error_rate = 0.03;

        // Prior probabilities
        let prior_homref = self.config.mixture_fractions[0].ln();
        let prior_het = self.config.mixture_fractions[1].ln();
        let prior_homalt = self.config.mixture_fractions[2].ln();

        let total = ref_reads + allele1_reads + allele2_reads;
        if total == 0 {
            return [0.0, 0.0, 0.0];
        }

        // 0/0: Both haplotypes are REF
        // Expect: mostly ref_reads, few allele1/allele2 (from errors)
        let ll_00 = {
            let binom = Binomial::new(1.0 - 2.0 * error_rate, total).unwrap();
            binom.ln_pmf(ref_reads)
        };

        // 0/1: One haplotype REF, one haplotype carries allele1
        // Possibility A: hap1=REF, hap2=allele1
        // Expect: ~50% ref_reads, ~50% allele1_reads, ~0% allele2_reads
        let ll_01_a = {
            // Model as trinomial, but use sequential binomials
            // First: P(allele2_reads | should be ~0)
            let p_error_allele2 = error_rate;
            let ll_allele2 = if total > 0 {
                Binomial::new(p_error_allele2, total)
                    .unwrap()
                    .ln_pmf(allele2_reads)
            } else {
                0.0
            };

            // Second: P(allele1_reads | remaining reads should split ~50/50 with ref)
            let remaining = ref_reads + allele1_reads;
            let ll_allele1 = if remaining > 0 {
                Binomial::new(0.5, remaining).unwrap().ln_pmf(allele1_reads)
            } else {
                0.0
            };

            ll_allele2 + ll_allele1
        };

        // Possibility B: hap1=REF, hap2=allele2
        //   Expect: ~50% ref_reads, ~0% allele1_reads, ~50% allele2_reads
        let ll_02_b = {
            let p_error_allele1 = error_rate;
            let ll_allele1 = if total > 0 {
                Binomial::new(p_error_allele1, total)
                    .unwrap()
                    .ln_pmf(allele1_reads)
            } else {
                0.0
            };

            let remaining = ref_reads + allele2_reads;
            let ll_allele2 = if remaining > 0 {
                Binomial::new(0.5, remaining).unwrap().ln_pmf(allele2_reads)
            } else {
                0.0
            };

            ll_allele1 + ll_allele2
        };

        // Marginalize over which allele is on which haplotype
        let max_het = ll_01_a.max(ll_02_b);
        let ll_01 = max_het + ((ll_01_a - max_het).exp() + (ll_02_b - max_het).exp()).ln();

        // 1/1: Both haplotypes carry alt alleles
        // Expect: mostly alt reads (allele1 + allele2), few ref_reads
        let ll_11 = {
            let alt_reads = allele1_reads + allele2_reads;
            let binom = Binomial::new(1.0 - error_rate, total).unwrap();
            binom.ln_pmf(alt_reads)
        };

        // Posterior = Prior + Likelihood
        [
            prior_homref + ll_00,
            prior_het + ll_01,
            prior_homalt + ll_11,
        ]
    }

    /// Calculates genotype quality (GQ) and sample quality (SQ)
    ///
    /// # Parameters
    /// - `gt_lplist`: The array from genotype_scores
    ///
    /// # Returns
    /// A tuple containing two floating-point values:
    /// - The first value is the genotype quality (GQ).
    /// - The second value is the sample quality (SQ).
    fn genotype_quals(&self, mut gt_lplist: [f64; 3]) -> (f64, f64) {
        // Convert from ln to log10
        gt_lplist
            .iter_mut()
            .for_each(|gt_lp| *gt_lp /= 10.0_f64.ln());

        // compute log10(sum of probabilities)
        let max_lp = gt_lplist.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let sum_exp: f64 = gt_lplist.iter().map(|&lp| 10.0_f64.powf(lp - max_lp)).sum();
        let log10_sum_exp = max_lp + sum_exp.log10();

        // Normalized log10 probabilities (capped at 0.0)
        let mut norm_log10_probs: Vec<f64> = gt_lplist
            .iter()
            .map(|&lp| f64::min(lp - log10_sum_exp, 0.0))
            .collect();

        let sq = f64::min((-10.0 * norm_log10_probs[0]).abs(), 1000.0);

        norm_log10_probs.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let gq = f64::min(-10.0 * (norm_log10_probs[1] - norm_log10_probs[0]), 1000.0);

        (self.calibrate_gq(gq), sq)
    }

    fn calibrate_gq(&self, raw_gq: f64) -> f64 {
        if self.config.calibration_table.is_empty() {
            return raw_gq;
        }

        // Linear interpolation
        if raw_gq <= self.config.calibration_table[0].0 {
            return self.config.calibration_table[0].1;
        }

        for i in 0..self.config.calibration_table.len() - 1 {
            let (x0, y0) = self.config.calibration_table[i];
            let (x1, y1) = self.config.calibration_table[i + 1];

            if raw_gq >= x0 && raw_gq <= x1 {
                // Linear interpolation
                let t = (raw_gq - x0) / (x1 - x0);
                return y0 + t * (y1 - y0);
            }
        }

        // Beyond table range
        self.config.calibration_table.last().unwrap().1
    }
}

impl Default for Genotyper {
    fn default() -> Self {
        Self {
            config: GenotyperConfig::default(),
        }
    }
}

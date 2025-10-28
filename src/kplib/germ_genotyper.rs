use ordered_float::OrderedFloat;
use statrs::{
    distribution::{Binomial, Discrete},
    function::gamma::ln_gamma,
};

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

/// Determines the genotype state based on coverage values for two alternate alleles.
/// The genotype state can be one of three: reference (Ref), heterozygous (Het), or homozygous (Hom).
/// If both coverage values are zero, the state is `Non`.
///
/// # Parameters
/// - `ref_cov`: The coverage value for the reference allele.
/// - `alt_cov`: The coverage value for the alternate allele.
///
/// # Returns
/// A `GTstate` enum value representing the genotype state
///
/// # Panics
/// This function will panic if an invalid state is encountered, which should be impossible under normal circumstances.
pub fn genotyper(ref_cov: u64, alt_cov: u64) -> GenotypeResult {
    let tot_cov = ref_cov + alt_cov;
    if tot_cov == 0 {
        return GenotypeResult {
            state: GTstate::Non,
            gq: 0.0,
            sq: 0.0,
        };
    }
    let scores = bino_genotype_scores(ref_cov, alt_cov);
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
    let (gq, sq) = genotype_quals(scores);
    GenotypeResult { state, gq, sq }
}

/// Calculates genotype scores for three possible genotypes (reference, heterozygous, homozygous)
/// based on the coverage values for two alternate alleles.
/// The scores are adjusted based on the total coverage to account for lower coverage scenarios.
///
/// # Parameters
/// - `ref_cov`: The coverage value for the reference allele.
/// - `alt_cov`: The coverage value for the alternate allele.
///
/// # Returns
/// An array of three floating-point values representing the log-probabilities for each genotype:
/// - The first value corresponds to the reference genotype.
/// - The second value corresponds to the heterozygous genotype.
/// - The third value corresponds to the homozygous genotype.
fn bino_genotype_scores(ref_cov: u64, alt_cov: u64) -> [f64; 3] {
    let error_rate = 0.01;

    // Prior probabilities (in log space)
    let prior_homref = 0.001_f64.ln();
    let prior_het = 0.75_f64.ln();
    let prior_homalt = 0.249_f64.ln();

    let n = ref_cov + alt_cov;
    if n == 0 {
        return [0.0, 0.0, 0.0];
    }

    // Create binomial distributions for each genotype
    let binom_homref = Binomial::new(error_rate, n).unwrap();
    let binom_het = Binomial::new(0.5, n).unwrap();
    let binom_homalt = Binomial::new(1.0 - error_rate, n).unwrap();

    // Calculate log-likelihoods
    let ll_homref = binom_homref.ln_pmf(alt_cov);
    let ll_het = binom_het.ln_pmf(alt_cov);
    let ll_homalt = binom_homalt.ln_pmf(alt_cov);

    // Posterior = Prior + Likelihood (in log space)
    [
        prior_homref + ll_homref,
        prior_het + ll_het,
        prior_homalt + ll_homalt,
    ]
}

fn __beta_binomial_ln_pmf(k: u64, n: u64, alpha: f64, beta: f64) -> f64 {
    let k = k as f64;
    let n = n as f64;

    // log C(n, k) = log(n!) - log(k!) - log((n-k)!)
    let log_binom_coef = ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0);

    // log B(k+α, n-k+β) = log Γ(k+α) + log Γ(n-k+β) - log Γ(n+α+β)
    let log_beta_num = ln_gamma(k + alpha) + ln_gamma(n - k + beta) - ln_gamma(n + alpha + beta);

    // log B(α, β) = log Γ(α) + log Γ(β) - log Γ(α+β)
    let log_beta_denom = ln_gamma(alpha) + ln_gamma(beta) - ln_gamma(alpha + beta);

    log_binom_coef + log_beta_num - log_beta_denom
}

fn __beta_genotype_scores(ref_cov: u64, alt_cov: u64) -> [f64; 3] {
    // IMPL OF BETA-BINOMIAL MODEL, TIES OUT WITH PYRO BETA-BINOMIAL IMPL WHEN HYPERPARAMETERS ARE SET CLOSE TO PYRO-FIT VALUES
    // TODO expose these as CLI parameters
    // for now, roughly set to typical values seen in HPRC samples fit with pyro implementation
    // // Hardy-Weinberg equilibrium:
    // let af = 0.3;
    // let p_homref = (1.0 - af).powi(2);  // 0.49 = 49% of sites
    // let p_het = 2.0 * af * (1.0 - af);   // 0.42 = 42% of sites
    // let p_homalt = af.powi(2);           // 0.09 = 9% of sites
    //let frac: &[f64] = &[0.2, 0.6, 0.2]; // mixture weights

    let total = ref_cov + alt_cov;

    if total == 0 {
        return [0.0, 0.0, 0.0]; // keep flat prior for missing GT to keep previous GQ < 5 threshold for LOWGQ filter
    }

    let frac: &[f64] = &[0.001, 0.75, 0.249]; // mixture weights
    let mu = &[0.005, 0.49, 0.99]; // beta-binomial means
    let nu = &[100.0, 46.90, 7.25]; // beta-binomial precisions

    // let coverage_factor = (total as f64 / 5.0).min(1.0);
    // let nu: Vec<f64> = nu_base.iter()
    // .map(|n| n * coverage_factor.max(0.2) * 0.50) // Minimum 20% of base nu
    // .collect();

    let alpha: Vec<f64> = mu.iter().zip(nu.iter()).map(|(m, n)| m * n).collect();
    let beta: Vec<f64> = mu
        .iter()
        .zip(nu.iter())
        .map(|(m, n)| (1.0 - m) * n)
        .collect();

    [
        frac[0].ln() + __beta_binomial_ln_pmf(alt_cov, total, alpha[0], beta[0]),
        frac[1].ln() + __beta_binomial_ln_pmf(alt_cov, total, alpha[1], beta[1]),
        frac[2].ln() + __beta_binomial_ln_pmf(alt_cov, total, alpha[2], beta[2]),
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
fn genotype_quals(mut gt_lplist: [f64; 3]) -> (f64, f64) {
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

    // Convert to linear probabilities for finding best genotype
    //let linear_probs: Vec<f64> = norm_log10_probs
    //.iter()
    //.map(|&log_p| 10.0_f64.powf(log_p))
    //.collect();

    // QUAL: Phred-scaled probability of homref
    // Since GL_P[0] = 10^(norm_log10_probs[0]), log10(GL_P[0]) = norm_log10_probs[0]
    let sq = f64::min((-10.0 * norm_log10_probs[0]).abs(), 1000.0);

    // GQ: Phred-scaled difference between best and second-best
    //let mut sorted = norm_log10_probs;
    norm_log10_probs.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let gq = f64::min(-10.0 * (norm_log10_probs[1] - norm_log10_probs[0]), 1000.0) / 10.0;

    // GQ: max of "probability it's NOT each genotype"
    //let gq_if_homref = -10.0 * (linear_probs[1] + linear_probs[2]).log10();
    //let gq_if_het = -10.0 * (linear_probs[0] + linear_probs[2]).log10();
    //let gq_if_homalt = -10.0 * (linear_probs[0] + linear_probs[1]).log10();

    // GQ: Use log-sum-exp to avoid underflow
    // For each genotype, compute log10(sum of other two probabilities)
    //let log_sum_not_homref = {
    //let m = norm_log10_probs[1].max(norm_log10_probs[2]);
    //let s = 10.0_f64.powf(norm_log10_probs[1] - m) + 10.0_f64.powf(norm_log10_probs[2] - m);
    //m + s.log10()
    //};

    //let log_sum_not_het = {
    //let m = norm_log10_probs[0].max(norm_log10_probs[2]);
    //let s = 10.0_f64.powf(norm_log10_probs[0] - m) + 10.0_f64.powf(norm_log10_probs[2] - m);
    //m + s.log10()
    //};

    //let log_sum_not_homalt = {
    //let m = norm_log10_probs[0].max(norm_log10_probs[1]);
    //let s = 10.0_f64.powf(norm_log10_probs[0] - m) + 10.0_f64.powf(norm_log10_probs[1] - m);
    //m + s.log10()
    //};

    //let gq_if_homref = -10.0 * log_sum_not_homref;
    //let gq_if_het = -10.0 * log_sum_not_het;
    //let gq_if_homalt = -10.0 * log_sum_not_homalt;

    //let gq = f64::min(gq_if_homref.max(gq_if_het).max(gq_if_homalt), 1000.0) / 10.0;

    (gq, sq)
}

fn __genotype_quals(mut gt_lplist: [f64; 3]) -> (f64, f64) {
    // Convert from ln to log10
    gt_lplist
        .iter_mut()
        .for_each(|gt_lp| *gt_lp /= 10.0_f64.ln());

    // Numerical stability: subtract max before exponentiating
    let max_lp = gt_lplist.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let probs: Vec<f64> = gt_lplist
        .iter()
        .map(|&lp| 10.0_f64.powf(lp - max_lp))
        .collect();
    let total: f64 = probs.iter().sum();

    let norm_probs: Vec<f64> = probs.iter().map(|p| p / total).collect();
    // SQ: Phred-scaled quality that it's not homref
    let sq = f64::min(-10.0 * (probs[0] / total).log10(), 100.0);

    // GQ: Phred-scaled difference between best and second-best
    // let mut sorted = gt_lplist;
    // sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
    // let delta = sorted[0] - sorted[1]; // Already in log10 space
    // let gq = f64::min(10.0 * delta, 100.0); // Just multiply by 10 for Phred

    // GQ: Maximum of "probability it's NOT each genotype"
    // For each genotype, sum the probabilities of the other two
    let gq_if_homref = -10.0 * (norm_probs[1] + norm_probs[2]).max(1e-10).log10();
    let gq_if_het = -10.0 * (norm_probs[0] + norm_probs[2]).max(1e-10).log10();
    let gq_if_homalt = -10.0 * (norm_probs[0] + norm_probs[1]).max(1e-10).log10();

    let gq = f64::max(gq_if_homref.max(gq_if_het).max(gq_if_homalt), 100.0);

    (gq, sq)
}

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
    let scores = genotype_scores(ref_cov, alt_cov);
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
    let (gq, sq) = genotype_quals(scores, tot_cov);
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
fn genotype_scores(ref_cov: u64, alt_cov: u64) -> [f64; 3] {
    let error_rate = 0.10;

    // Prior probabilities (in log space)
    let prior_homref = 0.001_f64.ln();
    let prior_het = 0.75_f64.ln();
    let prior_homalt = 0.249_f64.ln();

    let n = ref_cov + alt_cov;
    if n == 0 {
        return [1.0, 1.0, 1.0];
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

fn beta_binomial_ln_pmf(k: u64, n: u64, alpha: f64, beta: f64) -> f64 {
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

fn __betabinom_genotype_scores(ref_cov: u64, alt_cov: u64) -> [f64; 3] {
    // IMPL OF BETA-BINOMIAL MODEL, TIES OUT WITH PYRO BETA-BINOMIAL IMPL WHEN HYPERPARAMETERS ARE SET CLOSE TO PYRO-FIT VALUES
    // TODO expose these as CLI parameters
    // for now, roughly set to typical values seen in HPRC samples fit with pyro implementation
    // // Hardy-Weinberg equilibrium:
    // let af = 0.3;
    // let p_homref = (1.0 - af).powi(2);  // 0.49 = 49% of sites
    // let p_het = 2.0 * af * (1.0 - af);   // 0.42 = 42% of sites
    // let p_homalt = af.powi(2);           // 0.09 = 9% of sites
    //let frac: &[f64] = &[0.2, 0.6, 0.2]; // mixture weights

    let frac: &[f64] = &[0.001, 0.75, 0.249]; // mixture weights
    let mu: &[f64] = &[0.005, 0.49, 0.99]; // beta-binomial means
    let nu: &[f64] = &[30.0, 81.76, 11.74]; // beta-binomial precisions
    let alpha: Vec<f64> = mu.iter().zip(nu.iter()).map(|(m, n)| m * n).collect();
    let beta: Vec<f64> = mu
        .iter()
        .zip(nu.iter())
        .map(|(m, n)| (1.0 - m) * n)
        .collect();

    let total = ref_cov + alt_cov;
    if total == 0 {
        return [1.0, 1.0, 1.0]; // keep flat prior for missing GT to keep previous GQ < 5 threshold for LOWGQ filter
    }

    [
        frac[0].ln() + beta_binomial_ln_pmf(alt_cov, total, alpha[0], beta[0]),
        frac[1].ln() + beta_binomial_ln_pmf(alt_cov, total, alpha[1], beta[1]),
        frac[2].ln() + beta_binomial_ln_pmf(alt_cov, total, alpha[2], beta[2]),
    ]
}

/// Calculates genotype quality (GQ) and sample quality (SQ) based on the coverage values for reference and alternate alleles.
///
/// # Parameters
/// - `gt_lplist`: The array from genotype_scores
///
/// # Returns
/// A tuple containing two floating-point values:
/// - The first value is the genotype quality (GQ).
/// - The second value is the sample quality (SQ).
fn genotype_quals(mut gt_lplist: [f64; 3], total_cov: u64) -> (f64, f64) {
    // Convert from ln to log10
    gt_lplist
        .iter_mut()
        .for_each(|gt_lp| *gt_lp /= 10.0_f64.ln());

    // Convert to linear probabilities
    let probs: Vec<f64> = gt_lplist.iter().map(|&lp| 10.0_f64.powf(lp)).collect();
    let total: f64 = probs.iter().sum();

    // SQ: quality that it's not homref
    let sq = f64::min(-10.0 * (probs[0] / total).log10(), 100.0);

    // GQ: quality of best genotype call
    let best_prob = probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // Flat error rate for calibrating GQ based
    let prob_wrong = (total - best_prob) / total; // + 0.000001
    let gq = f64::min(-10.0 * prob_wrong.log10(), 100.0);

    (gq, sq)
}

fn __old_genotype_quals(mut gt_lplist: [f64; 3]) -> (f64, f64) {
    gt_lplist
        .iter_mut()
        .for_each(|gt_lp| *gt_lp /= 10.0_f64.ln());

    let mut gt_sum = 0.0;
    for gt in &gt_lplist {
        gt_sum += 10.0_f64.powf(*gt);
    }
    let gt_sum_log = gt_sum.log10();
    let sq = f64::min((-10.0 * (gt_lplist[0] - gt_sum_log)).abs(), 100.0);

    gt_lplist.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let best = gt_lplist[0];
    let second_best = gt_lplist[1];
    let gq = f64::min(-10.0 * (second_best - best), 100.0);

    (gq, sq)
}

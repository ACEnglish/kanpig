use ordered_float::OrderedFloat;
use rv::prelude::*;

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
/// - `alt1_cov`: The coverage value for the first alternate allele.
/// - `alt2_cov`: The coverage value for the second alternate allele.
///
/// # Returns
/// A `GTstate` enum value representing the genotype state
///
/// # Panics
/// This function will panic if an invalid state is encountered, which should be impossible under normal circumstances.
pub fn genotyper(alt1_cov: u64, alt2_cov: u64) -> GenotypeResult {
    if (alt1_cov + alt2_cov) == 0 {
        return GenotypeResult {
            state: GTstate::Non,
            gq: 0.0,
            sq: 0.0,
        };
    }
    let scores = genotype_scores(alt1_cov, alt2_cov);
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
/// - `alt1_cov`: The coverage value for the first alternate allele.
/// - `alt2_cov`: The coverage value for the second alternate allele.
///
/// # Returns
/// An array of three floating-point values representing the log-probabilities for each genotype:
/// - The first value corresponds to the reference genotype.
/// - The second value corresponds to the heterozygous genotype.
/// - The third value corresponds to the homozygous genotype.
fn genotype_scores(alt1_cov: u64, alt2_cov: u64) -> [f64; 3] {
    // IMPL OF BETA-BINOMIAL MODEL, TIES OUT WITH PYRO BETA-BINOMIAL IMPL WHEN HYPERPARAMETERS ARE SET CLOSE TO PYRO-FIT VALUES
    // TODO expose these as CLI parameters
    // for now, roughly set to typical values seen in HPRC samples fit with pyro implementation
    let frac: &[f64] = &[0.2, 0.6, 0.2]; // mixture weights
    let mu: &[f64] = &[0.005, 0.49, 0.99]; // beta-binomial means
    let nu: &[f64] = &[5.0, 50.0, 5.0]; // beta-binomial precisions
    let alpha: Vec<f64> = mu.iter().zip(nu.iter()).map(|(m, n)| m * n).collect();
    let beta: Vec<f64> = mu
        .iter()
        .zip(nu.iter())
        .map(|(m, n)| (1.0 - m) * n)
        .collect();

    let total = (alt1_cov + alt2_cov) as u32;
    if total == 0 {
        [1.0, 1.0, 1.0] // keep flat prior for missing GT to keep previous GQ < 5 threshold for LOWGQ filter
    } else {
        let beta_binom_homref = BetaBinomial::new(total, alpha[0], beta[0]).unwrap();
        let beta_binom_het = BetaBinomial::new(total, alpha[1], beta[1]).unwrap();
        let beta_binom_homalt = BetaBinomial::new(total, alpha[2], beta[2]).unwrap();

        [
            frac[0].ln() + beta_binom_homref.ln_pmf(&alt2_cov),
            frac[1].ln() + beta_binom_het.ln_pmf(&alt2_cov),
            frac[2].ln() + beta_binom_homalt.ln_pmf(&alt2_cov),
        ]
    }
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
fn genotype_quals(mut gt_lplist: [f64; 3]) -> (f64, f64) {
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

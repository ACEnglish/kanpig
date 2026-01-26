use std::collections::HashMap;

// Statistical functions
fn ln_gamma(x: f64) -> f64 {
    // Approximation of ln(Gamma(x)) using Stirling's formula for x > 1
    if x < 1.0 {
        return ln_gamma(x + 1.0) - x.ln();
    }
    (x - 0.5) * x.ln() - x + 0.5 * (2.0 * std::f64::consts::PI).ln()
}

fn ln_factorial(n: u32) -> f64 {
    if n <= 1 {
        0.0
    } else {
        ln_gamma(n as f64 + 1.0)
    }
}

fn ln_multinomial_pmf(counts: &[u32], probs: &[f64]) -> f64 {
    let n: u32 = counts.iter().sum();
    let mut log_prob = ln_factorial(n);

    for (count, prob) in counts.iter().zip(probs.iter()) {
        log_prob -= ln_factorial(*count);
        if *count > 0 && *prob > 0.0 {
            log_prob += (*count as f64) * prob.ln();
        } else if *count > 0 && *prob == 0.0 {
            return f64::NEG_INFINITY; // Impossible event
        }
    }
    log_prob
}

fn ln_beta_pdf(x: f64, alpha: f64, beta: f64) -> f64 {
    if x <= 0.0 || x >= 1.0 {
        return f64::NEG_INFINITY;
    }
    (alpha - 1.0) * x.ln() + (beta - 1.0) * (1.0 - x).ln() - ln_gamma(alpha) - ln_gamma(beta)
        + ln_gamma(alpha + beta)
}

#[derive(Debug, Clone)]
pub struct GenotypeHypothesis {
    pub germline_alleles: Vec<usize>, // indices of germline alleles (1 or 2)
    pub somatic_vafs: HashMap<usize, f64>, // VAFs for somatic alleleles
    pub observed_alleles: Vec<usize>, // indices of all alleles with coverage
}

#[derive(Debug, Clone)]
pub struct GenotypingResult {
    pub genotype: GenotypeHypothesis,
    pub log_posterior: f64,
    pub quality_score: f64,
    pub normalized_probabilities: Vec<f64>, // Probabilities of all hypotheses
}

/* Beta Distribution Basics
   The Beta distribution with parameters α (alpha) and β (beta) is defined on the interval [0, 1],
   making it perfect for modeling proportions like VAFs.
   How Alpha and Beta Shape the Distribution

   Mean: α / (α + β)
   Variance: αβ / [(α + β)²(α + β + 1)]
   Shape:

   α > 1, β > 1: Bell-shaped
   α < 1, β < 1: U-shaped
   α = β = 1: Uniform distribution

   somatic_vaf_prior_alpha: 1.0,    // α parameter
   somatic_vaf_prior_beta: 10.0,    // β parameter
   This gives you:

   Mean somatic VAF: 1.0 / (1.0 + 10.0) = 0.091 (~9%)
   Strong bias toward low VAFs: The distribution is heavily skewed toward 0
   Biological rationale: Most somatic mutations have low VAFs (5-20%), so this prior encodes that expectation

   Visual Interpretation
   With α=1, β=10:

   Peak probability near 0
   Rapidly decreasing as VAF increases
   Very low probability for VAFs > 0.3

   Tuning Guidelines
   More conservative (favor lower VAFs):
   alpha: 0.5, beta: 20.0  // Mean ~2.4%, very low VAF bias
   More permissive:
   alpha: 2.0, beta: 8.0   // Mean ~20%, allows higher somatic VAFs
   Uniform (no bias):
   alpha: 1.0, beta: 1.0   // Mean 50%, no preference
*/
pub struct MosaicGenotyper {
    pub error_rate: f64,
    pub somatic_vaf_prior_alpha: f64,
    pub somatic_vaf_prior_beta: f64,
    pub max_somatic_vaf: f64,
    pub min_depth_for_call: u32,
}

impl Default for MosaicGenotyper {
    fn default() -> Self {
        Self {
            error_rate: 0.001,
            somatic_vaf_prior_alpha: 1.0,
            somatic_vaf_prior_beta: 15.0,
            max_somatic_vaf: 0.2,
            min_depth_for_call: 1,
        }
    }
}

impl MosaicGenotyper {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_params(
        error_rate: f64,
        somatic_alpha: f64,
        somatic_beta: f64,
        max_somatic_vaf: f64,
        min_depth_for_call: u32,
    ) -> Self {
        Self {
            error_rate,
            somatic_vaf_prior_alpha: somatic_alpha,
            somatic_vaf_prior_beta: somatic_beta,
            max_somatic_vaf,
            min_depth_for_call,
        }
    }

    /// Generate all possible genotype hypotheses
    fn generate_hypotheses(&self, observed_alleles: &[usize]) -> Vec<GenotypeHypothesis> {
        let mut hypotheses = Vec::new();

        // Homozygous hypotheses
        for i in observed_alleles {
            let mut somatic_vafs = HashMap::new();
            for j in observed_alleles {
                if i != j {
                    somatic_vafs.insert(*j, 0.05);
                }
            }
            hypotheses.push(GenotypeHypothesis {
                germline_alleles: vec![*i],
                somatic_vafs,
                observed_alleles: Vec::new(),
            });
        }

        // Heterozygous hypotheses
        for (idx_i, &i) in observed_alleles.iter().enumerate() {
            for &j in &observed_alleles[idx_i + 1..] {
                let mut somatic_vafs = HashMap::new();
                for k in observed_alleles {
                    if *k != i && *k != j {
                        somatic_vafs.insert(*k, 0.05);
                    }
                }
                hypotheses.push(GenotypeHypothesis {
                    germline_alleles: vec![i, j],
                    somatic_vafs,
                    observed_alleles: Vec::new(),
                });
            }
        }

        hypotheses
    }

    /// Calculate expected VAFs without renormalizing
    fn calculate_expected_vafs(
        &self,
        hypothesis: &GenotypeHypothesis,
        num_alleles: usize,
    ) -> Vec<f64> {
        let num_germline = hypothesis.germline_alleles.len();

        // Calculate total somatic mass
        let somatic_mass: f64 = hypothesis.somatic_vafs.values().sum();

        // Remaining mass for germline alleles
        let germline_mass = (1.0 - somatic_mass).max(0.01); // Ensure positive

        let mut vafs = vec![0.0; num_alleles];

        // Distribute germline mass
        match num_germline {
            1 => {
                let idx = hypothesis.germline_alleles[0];
                vafs[idx] = germline_mass;
            }
            2 => {
                let share = germline_mass / 2.0;
                vafs[hypothesis.germline_alleles[0]] = share;
                vafs[hypothesis.germline_alleles[1]] = share;
            }
            _ => panic!("Invalid germline count"),
        }

        // Add somatic VAFs
        for (&idx, &vaf) in &hypothesis.somatic_vafs {
            vafs[idx] = vaf;
        }

        vafs
    }

    fn log_likelihood(&self, counts: &[u32], hypothesis: &GenotypeHypothesis) -> f64 {
        let expected_vafs = self.calculate_expected_vafs(hypothesis, counts.len());
        ln_multinomial_pmf(counts, &expected_vafs)
    }

    /// Improved prior using allele frequencies if available
    fn log_prior(&self, hypothesis: &GenotypeHypothesis, num_alleles: usize) -> f64 {
        let mut log_prior = 0.0;

        // Prior on somatic VAFs (Beta distribution)
        for &somatic_vaf in hypothesis.somatic_vafs.values() {
            if somatic_vaf > 0.0 && somatic_vaf < 1.0 {
                log_prior += ln_beta_pdf(
                    somatic_vaf,
                    self.somatic_vaf_prior_alpha,
                    self.somatic_vaf_prior_beta,
                );
            } else if somatic_vaf >= 1.0 {
                return f64::NEG_INFINITY;
            }
        }

        // Germline prior based
        log_prior += match hypothesis.germline_alleles.len() {
            1 => -(num_alleles as f64).ln(),
            2 => {
                let n_het = (num_alleles * (num_alleles - 1)) / 2;
                -(n_het as f64).ln()
            }
            _ => f64::NEG_INFINITY,
        };

        log_prior
    }

    /// Improved optimization with more granular search
    fn optimize_somatic_vafs(
        &self,
        counts: &[u32],
        hypothesis: &GenotypeHypothesis,
    ) -> (GenotypeHypothesis, f64) {
        let total_depth: u32 = counts.iter().sum();
        if total_depth == 0 {
            return (hypothesis.clone(), f64::NEG_INFINITY);
        }

        let somatic_indices: Vec<usize> = hypothesis.somatic_vafs.keys().copied().collect();

        if somatic_indices.is_empty() {
            let ll = self.log_likelihood(counts, hypothesis);
            return (hypothesis.clone(), ll);
        }

        let mut best_hypothesis = hypothesis.clone();
        let mut best_posterior = f64::NEG_INFINITY;

        let total_depth_f64 = total_depth as f64;

        // For each somatic allele, try multiple VAF candidates
        for allele_idx in &somatic_indices {
            let observed_vaf = counts[*allele_idx] as f64 / total_depth_f64;

            // Generate candidates: observed VAF + grid around it
            let mut candidates = vec![observed_vaf];

            // Add grid points
            for step in 0..=20 {
                let vaf = self.max_somatic_vaf * (step as f64 / 20.0);
                if vaf > 0.001 && vaf <= self.max_somatic_vaf {
                    candidates.push(vaf);
                }
            }

            // Try each candidate
            for &candidate_vaf in &candidates {
                let mut temp_hypothesis = best_hypothesis.clone();
                temp_hypothesis
                    .somatic_vafs
                    .insert(*allele_idx, candidate_vaf);

                let ll = self.log_likelihood(counts, &temp_hypothesis);
                let prior = self.log_prior(&temp_hypothesis, counts.len());
                let posterior = ll + prior;

                if posterior > best_posterior {
                    best_posterior = posterior;
                    best_hypothesis = temp_hypothesis;
                }
            }
        }

        (best_hypothesis, best_posterior)
    }

    /// Compute proper normalized quality score
    fn compute_quality_score(
        &self,
        best_idx: usize,
        all_log_posteriors: &[f64],
    ) -> (f64, Vec<f64>) {
        // Convert from ln to log10
        let all_log10: Vec<f64> = all_log_posteriors
            .iter()
            .map(|&lp| lp / 10.0_f64.ln())
            .collect();

        // Normalize using log-sum-exp
        let max_lp = all_log10.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let sum_exp: f64 = all_log10.iter().map(|&lp| 10.0_f64.powf(lp - max_lp)).sum();
        let log10_sum_exp = max_lp + sum_exp.log10();

        // Normalized log10 probabilities
        let norm_log10: Vec<f64> = all_log10
            .iter()
            .map(|&lp| f64::min(lp - log10_sum_exp, 0.0))
            .collect();

        // Convert to linear probabilities
        let linear_probs: Vec<f64> = norm_log10.iter().map(|&lp| 10.0_f64.powf(lp)).collect();

        // Calculate P(wrong) = sum of all probabilities except best
        let prob_wrong: f64 = linear_probs
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != best_idx)
            .map(|(_, &p)| p)
            .sum();

        // GQ = -10 * log10(P_wrong)
        let gq = if prob_wrong > 0.0 {
            f64::min(-10.0 * prob_wrong.log10(), 1000.0)
        } else {
            1000.0
        };

        (gq / 10.0, linear_probs)
    }

    pub fn genotype(&self, allele_counts: &[u32]) -> Option<GenotypingResult> {
        let total_depth: u32 = allele_counts.iter().sum();

        if total_depth < self.min_depth_for_call {
            return None;
        }

        let non_zero_alleles: Vec<usize> = allele_counts
            .iter()
            .enumerate()
            .filter(|(_, &count)| count > 0)
            .map(|(idx, _)| idx)
            .collect();

        if non_zero_alleles.is_empty() {
            return None;
        }

        let hypotheses = self.generate_hypotheses(&non_zero_alleles);
        let mut optimized_hypotheses = Vec::new();
        let mut all_log_posteriors = Vec::new();

        for hypothesis in hypotheses {
            // Optimize somatic VAFs
            let (optimized_hypothesis, log_likelihood) =
                self.optimize_somatic_vafs(allele_counts, &hypothesis);
            let log_prior = self.log_prior(&optimized_hypothesis, allele_counts.len());
            let log_posterior = log_likelihood + log_prior;

            optimized_hypotheses.push(optimized_hypothesis);
            all_log_posteriors.push(log_posterior);
        }

        // Find best hypothesis
        let best_idx = all_log_posteriors
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(idx, _)| idx)?;

        let mut best_genotype = optimized_hypotheses[best_idx].clone();
        let best_log_posterior = all_log_posteriors[best_idx];

        // Compute quality score with proper normalization
        let (quality_score, normalized_probs) =
            self.compute_quality_score(best_idx, &all_log_posteriors);

        // Filter to observed alleles
        best_genotype
            .germline_alleles
            .retain(|&x| non_zero_alleles.contains(&x));
        best_genotype.observed_alleles = non_zero_alleles.clone();
        best_genotype
            .somatic_vafs
            .retain(|&k, _| non_zero_alleles.contains(&k));

        Some(GenotypingResult {
            genotype: best_genotype,
            log_posterior: best_log_posterior,
            quality_score,
            normalized_probabilities: normalized_probs,
        })
    }
}

// Example usage and tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_het() {
        let genotyper = MosaicGenotyper::default();
        let counts = vec![50, 50, 0]; // Clear het between allele 0 and 1
        println!("Covs: {:?}", counts);
        if let Some(result) = genotyper.genotype(&counts) {
            println!("Genotype: {:?}", result.genotype.germline_alleles);
            println!("Quality Score: {:.2}", result.quality_score);
            println!("Normalized Probs: {:?}", result.normalized_probabilities);
            assert!(result.quality_score > 10.0); // Should be confident
        }
    }

    #[test]
    fn test_somatic() {
        let genotyper = MosaicGenotyper::default();
        let counts = vec![85, 10, 5]; // Mostly allele 0, some allele 1 (somatic?)
        println!("Covs: {:?}", counts);

        if let Some(result) = genotyper.genotype(&counts) {
            println!("Genotype: {:?}", result.genotype.germline_alleles);
            println!("Somatic VAFs: {:?}", result.genotype.somatic_vafs);
            println!("Quality Score: {:.2}", result.quality_score);
        }
    }

    #[test]
    fn test_somatic2() {
        let genotyper = MosaicGenotyper::default();
        let counts = vec![0, 34, 22, 10, 5]; // Mostly alt 0, some allele 1 (somatic?)
        println!("Covs: {:?}", counts);

        if let Some(result) = genotyper.genotype(&counts) {
            println!("Genotype: {:?}", result.genotype.germline_alleles);
            println!("Somatic VAFs: {:?}", result.genotype.somatic_vafs);
            println!("Quality Score: {:.2}", result.quality_score);
        }
    }
}

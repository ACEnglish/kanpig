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
    pub somatic_vafs: HashMap<usize, f64>, // VAFs for somatic alleles
    pub observed_alleles: Vec<usize>, // indices of all alleles with non-zero coverage
}

#[derive(Debug, Clone)]
pub struct GenotypingResult {
    pub genotype: GenotypeHypothesis,
    pub log_posterior: f64,
    pub quality_score: f64,
}

/*
* Beta Distribution Basics
   The Beta distribution with parameters α (alpha) and β (beta) is defined on the interval [0, 1], making it perfect for modeling proportions like VAFs.
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
    pub somatic_vaf_prior_alpha: f64, // Beta prior parameters for somatic VAFs
    pub somatic_vaf_prior_beta: f64,
    pub max_somatic_vaf: f64,
    pub min_depth_for_call: u32,
}

impl Default for MosaicGenotyper {
    fn default() -> Self {
        Self {
            error_rate: 0.001,            // don't know what these do
            somatic_vaf_prior_alpha: 1.0, // TODO: maybe a parameter
            somatic_vaf_prior_beta: 15.0, // Prior favoring low VAFs
            max_somatic_vaf: 0.2,         // TODO: Probably need this as a param
            min_depth_for_call: 1,        // TODO: I don't know if I need this as a param
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
    fn generate_hypotheses(&self, num_alleles: usize) -> Vec<GenotypeHypothesis> {
        let mut hypotheses = Vec::new();

        // Homozygous hypotheses (one germline allele)
        for i in 0..num_alleles {
            let mut somatic_vafs = HashMap::new();
            for j in 0..num_alleles {
                if i != j {
                    // Start with a reasonable somatic VAF estimate
                    somatic_vafs.insert(j, 0.05);
                }
            }
            hypotheses.push(GenotypeHypothesis {
                germline_alleles: vec![i],
                somatic_vafs,
                observed_alleles: Vec::new(), // Will be populated later
            });
        }

        // Heterozygous hypotheses (two germline alleles)
        for i in 0..num_alleles {
            for j in (i + 1)..num_alleles {
                let mut somatic_vafs = HashMap::new();
                for k in 0..num_alleles {
                    if k != i && k != j {
                        somatic_vafs.insert(k, 0.05);
                    }
                }
                hypotheses.push(GenotypeHypothesis {
                    germline_alleles: vec![i, j],
                    somatic_vafs,
                    observed_alleles: Vec::new(), // Will be populated later
                });
            }
        }

        hypotheses
    }

    /// Calculate expected VAFs under a given genotype hypothesis
    fn calculate_expected_vafs(
        &self,
        hypothesis: &GenotypeHypothesis,
        num_alleles: usize,
    ) -> Vec<f64> {
        let mut vafs = vec![self.error_rate; num_alleles];

        match hypothesis.germline_alleles.len() {
            1 => {
                // Homozygous
                let germline_idx = hypothesis.germline_alleles[0];
                vafs[germline_idx] = 1.0 - self.error_rate * (num_alleles - 1) as f64;
            }
            2 => {
                // Heterozygous
                let total_error = self.error_rate * (num_alleles - 2) as f64;
                let germline_vaf = (1.0 - total_error) / 2.0;
                vafs[hypothesis.germline_alleles[0]] = germline_vaf;
                vafs[hypothesis.germline_alleles[1]] = germline_vaf;
            }
            _ => panic!("Invalid number of germline alleles"),
        }

        // Set somatic VAFs
        for (&allele_idx, &somatic_vaf) in &hypothesis.somatic_vafs {
            vafs[allele_idx] = somatic_vaf;
        }

        // Normalize to ensure they sum to 1
        let sum: f64 = vafs.iter().sum();
        if sum > 0.0 {
            for vaf in &mut vafs {
                *vaf /= sum;
            }
        }

        vafs
    }

    /// Calculate log likelihood of observed data under hypothesis
    fn log_likelihood(&self, counts: &[u32], hypothesis: &GenotypeHypothesis) -> f64 {
        let expected_vafs = self.calculate_expected_vafs(hypothesis, counts.len());
        ln_multinomial_pmf(counts, &expected_vafs)
    }

    /// Calculate log prior probability of hypothesis
    fn log_prior(&self, hypothesis: &GenotypeHypothesis) -> f64 {
        let mut log_prior = 0.0;

        // Prior on somatic VAFs (Beta distribution)
        for &somatic_vaf in hypothesis.somatic_vafs.values() {
            log_prior += ln_beta_pdf(
                somatic_vaf,
                self.somatic_vaf_prior_alpha,
                self.somatic_vaf_prior_beta,
            );
        }

        // Simple uniform prior on germline genotypes for now
        // Could be improved with population allele frequencies
        log_prior += match hypothesis.germline_alleles.len() {
            1 => -1.0_f64.ln(), // Log uniform over homozygous states
            2 => -2.0_f64.ln(), // Log uniform over heterozygous states
            _ => f64::NEG_INFINITY,
        };

        log_prior
    }

    /// Optimize somatic VAFs for a given hypothesis using simple grid search
    fn optimize_somatic_vafs(
        &self,
        counts: &[u32],
        hypothesis: &GenotypeHypothesis,
    ) -> (GenotypeHypothesis, f64) {
        let total_depth: u32 = counts.iter().sum();
        if total_depth == 0 {
            return (hypothesis.clone(), f64::NEG_INFINITY);
        }

        let mut optimized_hypothesis = hypothesis.clone();

        // Simple optimization: try different VAF values for each somatic allele
        let vaf_candidates: Vec<f64> = (1..=20).map(|i| (i as f64) * 0.01).collect(); // 0.01 to 0.20

        // Get all somatic allele indices first to avoid borrowing issues
        let somatic_allele_indices: Vec<usize> = hypothesis.somatic_vafs.keys().copied().collect();

        for allele_idx in somatic_allele_indices {
            let observed_vaf = counts[allele_idx] as f64 / total_depth as f64;
            let current_vaf = optimized_hypothesis.somatic_vafs[&allele_idx];
            let mut best_vaf = current_vaf;
            let mut best_local_likelihood = f64::NEG_INFINITY;

            // Try different VAF values
            for &candidate_vaf in &vaf_candidates {
                if candidate_vaf > self.max_somatic_vaf {
                    break;
                }

                // Create temporary hypothesis with this VAF
                let mut temp_hypothesis = optimized_hypothesis.clone();
                temp_hypothesis
                    .somatic_vafs
                    .insert(allele_idx, candidate_vaf);
                let likelihood = self.log_likelihood(counts, &temp_hypothesis);

                if likelihood > best_local_likelihood {
                    best_local_likelihood = likelihood;
                    best_vaf = candidate_vaf;
                }
            }

            // Also try the observed VAF if it's reasonable
            if observed_vaf > 0.0 && observed_vaf <= self.max_somatic_vaf {
                let mut temp_hypothesis = optimized_hypothesis.clone();
                temp_hypothesis
                    .somatic_vafs
                    .insert(allele_idx, observed_vaf);
                let likelihood = self.log_likelihood(counts, &temp_hypothesis);
                if likelihood > best_local_likelihood {
                    best_vaf = observed_vaf;
                }
            }

            // Update the optimized hypothesis with the best VAF
            optimized_hypothesis
                .somatic_vafs
                .insert(allele_idx, best_vaf);
        }

        let final_likelihood = self.log_likelihood(counts, &optimized_hypothesis);
        (optimized_hypothesis, final_likelihood)
    }

    /// Main genotyping function
    pub fn genotype(&self, allele_counts: &[u32]) -> Option<GenotypingResult> {
        let total_depth: u32 = allele_counts.iter().sum();

        if total_depth < self.min_depth_for_call {
            return None;
        }

        // Filter out alleles with zero counts for hypothesis generation
        let non_zero_alleles: Vec<usize> = allele_counts
            .iter()
            .enumerate()
            .filter(|(_, &count)| count > 0)
            .map(|(idx, _)| idx)
            .collect();
        // Shouldn't happen
        if non_zero_alleles.is_empty() {
            return None;
        }

        let hypotheses = self.generate_hypotheses(allele_counts.len());
        let mut best_hypothesis = None;
        let mut best_log_posterior = f64::NEG_INFINITY;
        let mut second_best_log_posterior = f64::NEG_INFINITY;

        for hypothesis in hypotheses {
            // Skip hypotheses where germline alleles have zero coverage
            if hypothesis
                .germline_alleles
                .iter()
                .any(|&idx| allele_counts[idx] == 0)
            {
                continue;
            }

            // Optimize somatic VAFs
            let (optimized_hypothesis, log_likelihood) =
                self.optimize_somatic_vafs(allele_counts, &hypothesis);
            let log_prior = self.log_prior(&optimized_hypothesis);
            let log_posterior = log_likelihood + log_prior;

            if log_posterior > best_log_posterior {
                second_best_log_posterior = best_log_posterior;
                best_log_posterior = log_posterior;
                best_hypothesis = Some(optimized_hypothesis);
            } else if log_posterior > second_best_log_posterior {
                second_best_log_posterior = log_posterior;
            }
        }

        best_hypothesis.map(|mut genotype| {
            // Filter germline alleles to only include observed ones
            genotype
                .germline_alleles
                .retain(|&x| non_zero_alleles.contains(&x));

            // Set observed alleles
            genotype.observed_alleles = non_zero_alleles.clone();

            // Filter somatic VAFs to only include observed alleles
            genotype
                .somatic_vafs
                .retain(|&k, _| non_zero_alleles.contains(&k));

            // Calculate quality score (difference in log posterior)
            let quality_score = if second_best_log_posterior.is_finite() {
                (best_log_posterior - second_best_log_posterior) / (10.0_f64.ln() / 10.0)
            // Convert to Phred-like scale
            } else {
                100.0 // Very high confidence if only one viable hypothesis
            };

            GenotypingResult {
                genotype,
                log_posterior: best_log_posterior,
                quality_score,
            }
        })
    }
}

// Example usage and testing
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_heterozygous() {
        let genotyper = MosaicGenotyper::new();

        // Simulate het with two alleles at ~50% each
        let counts = vec![45, 55, 2, 1]; // Two main alleles + noise

        if let Some(result) = genotyper.genotype(&counts) {
            println!("Genotype result: {:#?}", result);
            assert_eq!(result.genotype.germline_alleles.len(), 2);
            assert!(result.genotype.germline_alleles.contains(&0));
            assert!(result.genotype.germline_alleles.contains(&1));
        } else {
            panic!("Should have produced a genotype call");
        }
    }

    #[test]
    fn test_homozygous_with_somatic() {
        let genotyper = MosaicGenotyper::new();

        // Simulate homozygous ref with somatic variants
        let counts = vec![90, 0, 8, 3]; // Dominant allele + somatic variants

        if let Some(result) = genotyper.genotype(&counts) {
            println!("Genotype result: {:#?}", result);
            assert_eq!(result.genotype.germline_alleles.len(), 1);
            assert_eq!(result.genotype.germline_alleles[0], 0);
            assert!(!result.genotype.somatic_vafs.is_empty());
        }
    }

    #[test]
    fn test_just_germline_het() {
        let mut genotyper = MosaicGenotyper::new();
        genotyper.somatic_vaf_prior_beta = 0.15;

        let counts = vec![47, 22, 4]; // Het with a tiny bit of noise

        if let Some(result) = genotyper.genotype(&counts) {
            println!("Genotype result: {:#?}", result);
            assert_eq!(result.genotype.germline_alleles.len(), 2);
            assert_eq!(result.genotype.germline_alleles[0], 1);
            assert!(result.genotype.somatic_vafs.is_empty());
        }
    }
}

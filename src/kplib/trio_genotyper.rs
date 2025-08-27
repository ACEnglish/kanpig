use ndarray::{Array1, Array2, Axis};
use std::f64;

/// Log-likelihood of observing k successes in n trials with probability p.
/// Simplified Bernoulli/Binomial form (no factorial terms).
fn log_binomial(k: usize, n: usize, p: f64) -> f64 {
    if p <= 0.0 || p >= 1.0 {
        if (k > 0 && p == 0.0) || (k < n && p == 1.0) {
            f64::NEG_INFINITY
        } else {
            0.0
        }
    } else {
        (k as f64) * p.ln() + ((n - k) as f64) * (1.0 - p).ln()
    }
}

/// Compute log-likelihood for a sample given haplotypes h1, h2.
fn compute_log_likelihood(
    h1: usize,
    h2: usize,
    sample: usize,
    n_clusters: usize,
    total_reads: usize,
    cluster_data: &Array2<usize>,
    qual: &[f64],
) -> f64 {
    if total_reads == 0 {
        return f64::NEG_INFINITY;
    }

    // expected distribution over clusters
    let mut expected: Vec<f64> = Vec::with_capacity(n_clusters);
    let mut total_exp = 0.0;
    for c in 0..n_clusters {
        let val = if c == h1 && c == h2 {
            1.0
        } else if c == h1 || c == h2 {
            0.5
        } else {
            1e-3
        };
        total_exp += val;
        expected.push(val);
    }

    // normalize
    for val in expected.iter_mut() {
        *val /= total_exp;
    }

    // compute likelihood
    let mut log_lik = 0.0;
    for c in 0..n_clusters {
        let obs = cluster_data[[c, sample]];
        let exp = expected[c];
        let q = qual[c];
        log_lik += log_binomial(obs, total_reads, exp) + (q + 1e-3f64).ln();
    }

    log_lik
}

/// Mendelian prior: probability that child genotype is consistent with parent genotypes.
fn mendelian_prior(gt_p: (usize, usize), gt_f: (usize, usize), gt_m: (usize, usize)) -> f64 {
    // get unique alleles (parents are diploid, so max 2 alleles each)
    let f_alleles = [gt_f.0, gt_f.1];
    let m_alleles = [gt_m.0, gt_m.1];

    // canonicalize child genotype
    let mut child = [gt_p.0, gt_p.1];
    child.sort_unstable();
    let gc_sorted = (child[0], child[1]);

    // build possible child genotypes
    let mut possible = Vec::new();
    for &a in &f_alleles {
        for &b in &m_alleles {
            let mut pair = [a, b];
            pair.sort_unstable();
            let tup = (pair[0], pair[1]);
            if !possible.contains(&tup) {
                possible.push(tup);
            }
        }
    }

    if possible.contains(&gc_sorted) {
        (1.0 / (possible.len() as f64)).ln()
    } else {
        (1e-6f64).ln()
    }
}

/// Compute phred-scaled genotype quality from posteriors.
/// Input: slice of (genotype, log-posterior), sorted descending by log-posterior.
fn __genotype_quality<T>(posteriors: &[(T, f64)]) -> i32 {
    if posteriors.is_empty() {
        return 0;
    }

    let best_logp = posteriors[0].1;
    if posteriors.len() == 1 {
        return 100;
    }

    let next_best_logp = posteriors[1].1;
    let delta = best_logp - next_best_logp;

    let prob_best = 1.0 / (1.0 + (-delta).exp());

    let gq = -10.0 * (1.0 - prob_best).max(1e-10).log10();

    gq.round().min(100.0) as i32
}

/// Joint trio genotyping: returns best genotype tuple (child, father, mother)
pub fn trio_genotyper(read_counts: &Array2<usize>, qual: &[f64]) -> [[usize; 2]; 3] {
    let n_clusters = read_counts.nrows();

    // generate all diploid genotypes
    let mut genotypes = Vec::new();
    for i in 0..n_clusters {
        for j in i..n_clusters {
            genotypes.push((i, j));
        }
    }

    let mut trio_posteriors: Vec<([[usize; 2]; 3], f64)> = Vec::new();

    let n_clusters = read_counts.nrows();
    let total_reads: Array1<usize> = read_counts.sum_axis(Axis(0));

    for &gt_f in &genotypes {
        for &gt_m in &genotypes {
            for &gt_p in &genotypes {
                let ll_c = compute_log_likelihood(
                    gt_p.0,
                    gt_p.1,
                    0,
                    n_clusters,
                    total_reads[0],
                    read_counts,
                    qual,
                );
                let ll_f = compute_log_likelihood(
                    gt_f.0,
                    gt_f.1,
                    1,
                    n_clusters,
                    total_reads[1],
                    read_counts,
                    qual,
                );
                let ll_m = compute_log_likelihood(
                    gt_m.0,
                    gt_m.1,
                    2,
                    n_clusters,
                    total_reads[2],
                    read_counts,
                    qual,
                );
                let lp = mendelian_prior(gt_p, gt_f, gt_m);
                let posterior = ll_f + ll_m + ll_c + lp;
                let gts_array = [[gt_p.0, gt_p.1], [gt_f.0, gt_f.1], [gt_m.0, gt_m.1]];
                trio_posteriors.push((gts_array, posterior));
            }
        }
    }

    let best = trio_posteriors
        .iter()
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .map(|(geno, _)| *geno)  // take the genotype part
        .unwrap();

    // let gq = genotype_quality(&trio_posteriors[..2]);
    best
}

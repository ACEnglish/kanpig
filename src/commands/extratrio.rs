use ndarray::Array2;
use std::collections::HashMap;
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
    cluster_data: &Array2<usize>,
    qual: &HashMap<usize, f64>,
) -> f64 {
    let n_clusters = cluster_data.nrows();

    // total reads in this sample
    let total_reads: usize = cluster_data
        .outer_iter()
        .map(|row| row[sample])
        .sum();

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
        // NOTE: qual was unused in your current Python version
        let _q = *qual.get(&c).unwrap_or(&0.0);
        log_lik += log_binomial(obs, total_reads, exp);
    }

    log_lik
}

/// Mendelian prior: probability that child genotype is consistent with parent genotypes.
fn mendelian_prior(Gf: (usize, usize), Gm: (usize, usize), Gc: (usize, usize)) -> f64 {
    // get unique alleles (parents are diploid, so max 2 alleles each)
    let f_alleles = [Gf.0, Gf.1];
    let m_alleles = [Gm.0, Gm.1];

    // canonicalize child genotype
    let mut child = [Gc.0, Gc.1];
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
fn genotype_quality<T>(posteriors: &[(T, f64)]) -> i32 {
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

/// Joint trio genotyping: returns best genotype tuple (child, father, mother) and quality.
fn trio_genotyper(
    read_counts: &Array2<usize>,
    qual: &HashMap<usize, f64>,
) -> (( (usize, usize), (usize, usize), (usize, usize) ), i32) {
    let n_clusters = read_counts.nrows();

    // generate all diploid genotypes
    let mut genotypes = Vec::new();
    for i in 0..n_clusters {
        for j in i..n_clusters {
            genotypes.push((i, j));
        }
    }

    // compute trio posteriors
    let mut trio_posteriors: Vec<(((usize, usize), (usize, usize), (usize, usize)), f64)> = Vec::new();

    for &Gf in &genotypes {
        for &Gm in &genotypes {
            for &Gc in &genotypes {
                let ll_m = compute_log_likelihood(Gm.0, Gm.1, 2, read_counts, qual);
                let ll_f = compute_log_likelihood(Gf.0, Gf.1, 1, read_counts, qual);
                let ll_c = compute_log_likelihood(Gc.0, Gc.1, 0, read_counts, qual);
                let lp = mendelian_prior(Gf, Gm, Gc);
                let posterior = ll_f + ll_m + ll_c + lp;
                trio_posteriors.push(((Gc, Gf, Gm), posterior));
            }
        }
    }

    // sort by posterior descending
    trio_posteriors.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // print top 5
    for i in trio_posteriors.iter().take(5) {
        println!("{:?}", i);
    }

    let best = &trio_posteriors[0];
    let gq = genotype_quality(&trio_posteriors[..2]);
    (best.0, gq)
}

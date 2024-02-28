use simsimd::SimSIMD;

/// Cosine similarity
pub fn cosinesim(a: &[f32], b: &[f32]) -> f32 {
    (1.0 - f32::cosine(a, b).unwrap()).abs()
}

/// Experimental jaccard-like similarity - it actually might just be a normalized euclidean
/// distance
pub fn seqsim(a: &[f32], b:&[f32]) -> f32 {
    println!("comparing {:?} {:?}", a, b);
    let deno = a.iter().zip(b.iter()).map(|(x,y)| x.abs() + y.abs()).sum::<f32>();
    if deno == 0.0 { // no kmers
        println!("nokmers");
        return 0.0;
    }

    let neum = a.iter().zip(b.iter()).map(|(x,y)| (x - y).abs()).sum::<f32>();
    if neum == 0.0 { // identical
        println!("identical");
        return 1.0;
    }
    let ret = 1.0 - (neum / deno);
    // Argument for 1.0 - 4446660
    // 4809 4810
    // path size 4810 sig p:1 t:1
    // 17 over 9613 equals 0.9982316

    // Argument against 1.0 - 6170297
    //-332 -332
    // path size -332 sig p:-1 t:-1
    // 659 over 659 equals 0
    println!("{} over {} equals {}", neum, deno, ret);
    ret
}

/// Weighted cosine similarity
/// Cosine similarity seems to be influenced by large outlier values. By weighing by frequency,
/// long sequences' cosine similarity behaves more consistently with small sequences
pub fn weighted_cosinesim(a: &[f32], b: &[f32]) -> f32 {
    let cos: f32 = cosinesim(a, b);
    let dist: f32 = (1.0 / a.len() as f32)
        * a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs())
            .sum::<f32>();

    cos / (dist.powi(2) + cos)
}

/// size similarity of two variant sizes
/// sizes must be positive
pub fn sizesim(size_a: u64, size_b: u64) -> f32 {
    if ((size_a == 0) || (size_b == 0)) && size_a == size_b {
        return 1.0;
    }
    std::cmp::max(std::cmp::min(size_a, size_b), 1) as f32
        / std::cmp::max(std::cmp::max(size_a, size_b), 1) as f32
}

/// do two intervals overlap
pub fn overlaps(s1: u64, e1: u64, s2: u64, e2: u64) -> bool {
    std::cmp::max(s1, s2) < std::cmp::min(e1, e2)
}

use std::f64::consts::LN_10;

/// Generate probabilities of genotypes
/// Smallest value is most likely genotype
/// I need this to return GT::REF,GT::HET,GT::HOM
/// Then, instead of relying on the hard threshold, we can check if genotyper(coverage, alt_cov) ==
/// GT::HOM
pub fn genotyper(tot_cov: f64, alt_cov: f64, priors: Option<&[f64]>) -> Option<Vec<f64>> {
    if tot_cov == 0.0 {
        return None;
    }

    let priors = priors.unwrap_or(&[0.05, 0.5, 0.95]);

    let mut gt_list = Vec::new();

    let total = tot_cov;
    let alt = alt_cov;
    let non_alt = total - alt;

    for &p_alt in priors {
        let mut comb = log_choose(total, alt);
        comb += alt * p_alt.ln();
        comb += non_alt * (1.0 - p_alt).ln();
        gt_list.push(comb);
    }

    Some(gt_list)
}

fn log_choose(n: f64, k: f64) -> f64 {
    let mut r = 0.0;
    let mut n = n;
    let mut k = k;

    if k * 2.0 > n {
        k = n - k;
    }

    for d in 1..=k as u64 {
        r += n.ln();
        r -= d as f64;
        n -= 1.0;
    }

    r / LN_10 // Convert to base 10 logarithm
}

use ordered_float::OrderedFloat;

/// Canberra distance of featurized kmers
pub fn seqsim(a: &[f32], b: &[f32], mink: f32) -> f32 {
    let mut deno: f32 = 0.0;
    let mut neum: f32 = 0.0;
    for (x, y) in a.iter().zip(b.iter()) {
        let d = x.abs() + y.abs();
        if d <= mink {
            continue;
        }
        deno += d;

        neum += (x - y).abs();
    }

    // no kmers
    if deno == 0.0 {
        return 0.0;
    }

    // identical
    if neum == 0.0 {
        return 1.0;
    }

    1.0 - (neum / deno)
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

#[derive(Debug, PartialEq, Eq)]
pub enum GTstate {
    Ref,
    Het,
    Hom,
    Non,
    //Hemi should be a thing
}

/// Generate genotypes given observed allele coverages
pub fn genotyper(alt1_cov: f64, alt2_cov: f64) -> GTstate {
    if (alt1_cov + alt2_cov) == 0.0 {
        return GTstate::Non;
    }
    match genotype_scores(alt1_cov, alt2_cov)
        .iter()
        .enumerate()
        .max_by_key(|&(_, &x)| OrderedFloat(x))
        .map(|(i, _)| i)
    {
        Some(0) => GTstate::Ref,
        Some(1) => GTstate::Het,
        Some(2) => GTstate::Hom,
        _ => panic!("not possible"),
    }
}

/// Probabilities of each genotype given the allele coveages
fn genotype_scores(alt1_cov: f64, alt2_cov: f64) -> Vec<f64> {
    let p_alt: &[f64] = &[1e-3, 0.5, 0.9];

    let total = alt1_cov + alt2_cov;
    let log_combo = log_choose(total, alt2_cov);

    let lp_homref = log_combo + alt2_cov * p_alt[0].log10() + alt1_cov * (1.0 - p_alt[0]).log10();
    let lp_het = log_combo + alt2_cov * p_alt[1].log10() + alt1_cov * (1.0 - p_alt[1]).log10();
    let lp_homalt = log_combo + alt2_cov * p_alt[2].log10() + alt1_cov * (1.0 - p_alt[2]).log10();

    vec![lp_homref, lp_het, lp_homalt]
}

/// Genotype quality: confidence in the assigned genotype
/// Sample quality: confidence that there is non-reference present
pub fn genotype_quals(ref_cov: f64, alt_cov: f64) -> (f64, f64) {
    let gt_lplist = genotype_scores(ref_cov, alt_cov);
    let mut sorted_gt_lplist: Vec<(usize, f64)> =
        gt_lplist.iter().enumerate().map(|(i, &e)| (i, e)).collect();
    sorted_gt_lplist.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    let best = sorted_gt_lplist[0];
    let second_best = sorted_gt_lplist[1];
    let gq = f64::min(-10.0 * (second_best.1 - best.1), 100.0);

    let mut gt_sum = 0.0;
    for gt in &gt_lplist {
        gt_sum += 10.0_f64.powf(*gt);
    }

    let gt_sum_log = gt_sum.log10();
    let sq = f64::min((-10.0 * (gt_lplist[0] - gt_sum_log)).abs(), 100.0);
    (gq, sq)
}

/// Helper function for genotype_scores
fn log_choose(n: f64, k: f64) -> f64 {
    let mut r = 0.0;
    let mut n = n;
    let mut k = k;

    if k * 2.0 > n {
        k = n - k;
    }

    for d in 1..=(k as i32) {
        r += n.log10();
        r -= d as f64;
        n -= 1.0;
    }

    r
}

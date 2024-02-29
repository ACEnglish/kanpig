use ordered_float::OrderedFloat;
use simsimd::SimSIMD;
use std::f64::consts::LN_10;

/// Cosine similarity
pub fn cosinesim(a: &[f32], b: &[f32]) -> f32 {
    (1.0 - f32::cosine(a, b).unwrap()).abs()
}

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


#[derive(Debug, PartialEq, Eq)]
pub enum GTstate {
    Ref,
    Het,
    Hom,
    Non,
    Unk,
    //Hemi should be a thing
}

/// Generate probabilities of genotypes
/// Smallest value is most likely genotype
pub fn genotyper(alt1_cov: f64, alt2_cov: f64) -> GTstate {
    let total = alt1_cov + alt2_cov;
    if total == 0.0 {
        return GTstate::Unk;
    }

    let priors: &[f64] = &[0.05, 0.5, 0.95];

    let mut gt_list = Vec::new();

    for &p_alt in priors {
        let mut comb = log_choose(total, alt1_cov);
        comb += alt2_cov * p_alt.ln();
        comb += alt1_cov * (1.0 - p_alt).ln();
        gt_list.push(comb);
    }

    // Some(gt_list) What I should be returning for a GQ
    let ret = match gt_list
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
    println!("GENOTYPER: {:?}", ret);
    ret
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

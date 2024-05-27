use ordered_float::OrderedFloat;

/// Computes the Canberra distance similarity between two featurized k-mer vectors.
/// The similarity is calculated as 1 minus the Canberra distance, providing a measure of similarity between 0 and 1.
///
/// # Parameters
/// - `a`: A slice of floating-point numbers representing the first k-mer vector.
/// - `b`: A slice of floating-point numbers representing the second k-mer vector.
/// - `mink`: A floating-point threshold below which differences are ignored.
///
/// # Returns
/// A floating-point value representing the similarity between the two vectors:
/// - 1.0 indicates identical vectors.
/// - 0.0 indicates no kmers or maximum dissimilarity.
pub fn seqsim(a: &[f32], b: &[f32], mink: f32) -> f32 {
    let mut deno: f32 = 0.0;
    let mut neum: f32 = 0.0;
    let mut total_d: f32;

    for (&x, &y) in a.iter().zip(b.iter()) {
        total_d = x.abs() + y.abs();
        if total_d > mink {
            deno += total_d;
            neum += (x - y).abs();
        }
    }

    if deno == 0.0 {
        return 0.0;
    }

    if neum == 0.0 {
        return 1.0;
    }

    1.0 - (neum / deno)
}

/// Computes size similarity
/// The similarity is defined as the ratio of the smaller size to the larger size,
/// with special handling for cases where either size is zero.
///
/// # Parameters
/// - `size_a`: The first size as a 64-bit unsigned integer.
/// - `size_b`: The second size as a 64-bit unsigned integer.
///
/// # Returns
/// A floating-point value representing the similarity score between the two sizes.
/// - If both sizes are zero, the function returns 1.0.
/// - Otherwise, the similarity is calculated as the ratio of the smaller size to the larger size.
pub fn sizesim(size_a: u64, size_b: u64) -> f32 {
    if size_a == size_b {
        return 1.0;
    }
    let min_size = size_a.min(size_b).max(1) as f32;
    let max_size = size_a.max(size_b).max(1) as f32;
    min_size / max_size
}

/// Determines if two intervals overlap.
/// Each interval is defined by a start and an end position.
/// The intervals overlap if the maximum of the start positions is less than the minimum of the end positions.
///
/// # Parameters
/// - `s1`: The start position of the first interval as a 64-bit unsigned integer.
/// - `e1`: The end position of the first interval as a 64-bit unsigned integer.
/// - `s2`: The start position of the second interval as a 64-bit unsigned integer.
/// - `e2`: The end position of the second interval as a 64-bit unsigned integer.
///
/// # Returns
/// A boolean value indicating whether the intervals overlap (`true`) or not (`false`).
pub fn overlaps(s1: u64, e1: u64, s2: u64, e2: u64) -> bool {
    std::cmp::max(s1, s2) < std::cmp::min(e1, e2)
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub enum GTstate {
    Ref,
    Het,
    Hom,
    Non,
    //Hemi should be a thing
}

/// Determines the genotype state based on coverage values for two alternate alleles.
/// The genotype state can be one of three: reference (Ref), heterozygous (Het), or homozygous (Hom).
/// If both coverage values are zero, the state is `Non`.
///
/// # Parameters
/// - `alt1_cov`: The coverage value for the first alternate allele as a floating-point number.
/// - `alt2_cov`: The coverage value for the second alternate allele as a floating-point number.
///
/// # Returns
/// A `GTstate` enum value representing the genotype state
///
/// # Panics
/// This function will panic if an invalid state is encountered, which should be impossible under normal circumstances.
pub fn genotyper(alt1_cov: f64, alt2_cov: f64) -> GTstate {
    if (alt1_cov + alt2_cov) == 0.0 {
        return GTstate::Non;
    }
    let ret = match genotype_scores(alt1_cov, alt2_cov)
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
    debug!("{} {} -> {:?}", alt1_cov, alt2_cov, ret);
    ret
}

/// Calculates genotype scores for three possible genotypes (reference, heterozygous, homozygous)
/// based on the coverage values for two alternate alleles.
/// The scores are adjusted based on the total coverage to account for lower coverage scenarios.
///
/// # Parameters
/// - `alt1_cov`: The coverage value for the first alternate allele as a floating-point number.
/// - `alt2_cov`: The coverage value for the second alternate allele as a floating-point number.
///
/// # Returns
/// An array of three floating-point values representing the log-probabilities for each genotype:
/// - The first value corresponds to the reference genotype.
/// - The second value corresponds to the heterozygous genotype.
/// - The third value corresponds to the homozygous genotype.
fn genotype_scores(alt1_cov: f64, alt2_cov: f64) -> [f64; 3] {
    // Needs to be more pure for lower coverage
    let p_alt: &[f64] = if alt1_cov + alt2_cov < 10.0 {
        &[1e-3, 0.55, 0.95]
    } else {
        &[1e-3, 0.50, 0.90]
    };

    let total = alt1_cov + alt2_cov;
    let log_combo = log_choose(total, alt2_cov);

    [
        log_combo + alt2_cov * p_alt[0].log10() + alt1_cov * (1.0 - p_alt[0]).log10(),
        log_combo + alt2_cov * p_alt[1].log10() + alt1_cov * (1.0 - p_alt[1]).log10(),
        log_combo + alt2_cov * p_alt[2].log10() + alt1_cov * (1.0 - p_alt[2]).log10(),
    ]
}

/// Calculates genotype quality (GQ) and sample quality (SQ) based on the coverage values for reference and alternate alleles.
///
/// # Parameters
/// - `ref_cov`: The coverage value for the reference allele as a floating-point number.
/// - `alt_cov`: The coverage value for the alternate allele as a floating-point number.
///
/// # Returns
/// A tuple containing two floating-point values:
/// - The first value is the genotype quality (GQ).
/// - The second value is the sample quality (SQ).
pub fn genotype_quals(ref_cov: f64, alt_cov: f64) -> (f64, f64) {
    let mut gt_lplist = genotype_scores(ref_cov, alt_cov);
    gt_lplist.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let best = gt_lplist[0];
    let second_best = gt_lplist[1];
    let gq = f64::min(-10.0 * (second_best - best), 100.0);

    let mut gt_sum = 0.0;
    for gt in &gt_lplist {
        gt_sum += 10.0_f64.powf(*gt);
    }

    let gt_sum_log = gt_sum.log10();
    let sq = f64::min((-10.0 * (gt_lplist[0] - gt_sum_log)).abs(), 100.0);
    (gq, sq)
}

/// Helper function for genotype_scores
const FACTORIAL_LIMIT: usize = 100;
lazy_static::lazy_static! {
    static ref LOG_FACTORIALS: Vec<f64> = {
        let mut log_factorials = vec![0.0; FACTORIAL_LIMIT + 1];
        let mut log_n_fact = 0.0;
        for (n, item) in log_factorials.iter_mut().enumerate().take(FACTORIAL_LIMIT + 1).skip(1) {
            log_n_fact += (n as f64).ln();
            *item = log_n_fact;
        }
        log_factorials
    };
}

fn log_choose(n: f64, k: f64) -> f64 {
    if n.is_infinite() || k.is_infinite() || n.is_nan() || k.is_nan() {
        return f64::NAN;
    }

    if k > n || k < 0.0 {
        return 0.0;
    }

    /*if k * 2.0 > n {
        k = n - k;
    }*/

    if n <= FACTORIAL_LIMIT as f64 {
        return LOG_FACTORIALS[n as usize]
            - LOG_FACTORIALS[k as usize]
            - LOG_FACTORIALS[(n - k) as usize];
    }

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

use crate::kplib::KmerVec;
/// Computes the Canberra distance similarity between two featurized k-mer vectors.
/// The similarity is calculated as 1 minus the Canberra distance, providing a measure of similarity between 0 and 1.
///
/// # Parameters
/// - `a`: A slice of floating-point numbers representing the first k-mer vector.
/// - `b`: A slice of floating-point numbers representing the second k-mer vector.
/// - `mink`: Threshold for minimum number of observations across vectors.
///
/// # Returns
/// A floating-point value representing the similarity between the two vectors:
/// - 1.0 indicates identical vectors.
/// - 0.0 indicates no kmers or maximum dissimilarity.
pub fn seqsim(a: &KmerVec, b: &KmerVec, mink: f32) -> f32 {
    let mut deno: f32 = 0.0;
    let mut neum: f32 = 0.0;

    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let (x, y) = match a[i].0.cmp(&b[j].0) {
            std::cmp::Ordering::Less => {
                i += 1;
                (a[i - 1].1, 0.0)
            }
            std::cmp::Ordering::Greater => {
                j += 1;
                (0.0, b[j - 1].1)
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
                (a[i - 1].1, b[j - 1].1)
            }
        };
        let total_d = x.abs() + y.abs();
        if total_d >= mink {
            deno += total_d;
            neum += (x - y).abs();
        }
    }

    // Unmatched tail entries — other side is 0.0
    for &(_, x) in &a[i..] {
        let total_d = x.abs();
        if total_d >= mink {
            deno += total_d;
            neum += total_d;
        }
    }
    for &(_, y) in &b[j..] {
        let total_d = y.abs();
        if total_d >= mink {
            deno += total_d;
            neum += total_d;
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

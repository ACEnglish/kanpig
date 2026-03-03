use crate::kplib::KmerVec;

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

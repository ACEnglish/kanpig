use simsimd::SimSIMD;

/// Cosine similarity
pub fn cosinesim(a: &[f32], b: &[f32]) -> f32 {
    1.0 - f32::cosine(a, b).unwrap()
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

use simsimd::SimSIMD;

/// Cosine similarity
pub fn cosinesim(a: &[f32], b: &[f32]) -> f32 {
    1.0 - f32::cosine(a, b).unwrap()
}

/// Weighted cosine similarity
pub fn weighted_cosinesim(a: &[f32], b: &[f32]) -> f32 {
    let cos: f32 = cosinesim(a, b);
    let dist: f32 = (1.0 / a.len() as f32)
        * a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).abs())
            .sum::<f32>();

    cos / (dist.powf(2.0) + cos)
}

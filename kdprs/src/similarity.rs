use simsimd::SimSIMD;

/// Cosine similarity
pub fn cosinesim(a: &Vec<f32>, b: &Vec<f32>) -> f32 {
    1.0 - f32::cosine(a, b).unwrap()
}

/// Weighted cosine similarity
pub fn weighted_cosinesim(a: &Vec<f32>, b: &Vec<f32>) -> f32 {
    let cos: f32 = cosinesim(a, b);
    let mut dist: f32 = a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).sum();
    dist = (1.0 / a.len() as f32) * dist;

    cos / (dist.powf(2.0) + cos)
}


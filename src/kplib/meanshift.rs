/// Mean Shift clustering algorithm for 1D data
pub struct MeanShift {
    pub bandwidth: Option<f64>,
    pub max_iter: usize,
    pub tol: f64,
    pub min_bin_freq: usize,
}

impl MeanShift {
    /// Create a new MeanShift instance
    pub fn new() -> Self {
        Self {
            bandwidth: None,
            max_iter: 300,
            tol: 1e-3,
            min_bin_freq: 1,
        }
    }

    /// Create a new MeanShift instance with specified bandwidth
    pub fn with_bandwidth(bandwidth: f64) -> Self {
        Self {
            bandwidth: Some(bandwidth),
            max_iter: 300,
            tol: 1e-3,
            min_bin_freq: 1,
        }
    }

    /// Set maximum iterations
    pub fn max_iter(mut self, max_iter: usize) -> Self {
        self.max_iter = max_iter;
        self
    }

    /// Set tolerance for convergence
    pub fn tolerance(mut self, tol: f64) -> Self {
        self.tol = tol;
        self
    }

    /// Set min_bin_freq
    pub fn min_size(mut self, minbin: usize) -> Self {
        self.min_bin_freq = minbin;
        self
    }

    /// Fit the Mean Shift algorithm to 1D data
    pub fn fit(&mut self, data: &[f64]) -> MeanShiftResult {
        if data.is_empty() {
            return MeanShiftResult {
                cluster_centers: vec![],
                labels: vec![],
                medoids: vec![],
            };
        }

        // Estimate bandwidth if not provided
        let bandwidth = match self.bandwidth {
            Some(bw) => bw,
            None => estimate_bandwidth(data, 0.3, None, 0).max(1.0),
        };

        // Use bin seeding like scikit-learn for better performance and consistency
        let seeds = self.get_bin_seeds(data, bandwidth, self.min_bin_freq);

        let mut centers = Vec::<f64>::new();

        // Apply mean shift to each seed
        for seed in seeds {
            let center = self.mean_shift_single_seed(data, seed, bandwidth);
            centers.push(center);
        }

        // Remove duplicate centers using bandwidth-based clustering
        // This is more similar to scikit-learn's approach
        centers.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut unique_centers = Vec::<f64>::new();

        for center in centers {
            let mut is_duplicate = false;
            for existing in &unique_centers {
                if (center - existing).abs() < bandwidth {
                    is_duplicate = true;
                    break;
                }
            }
            if !is_duplicate {
                unique_centers.push(center);
            }
        }

        // Assign labels to data points
        let labels = self.assign_labels(data, &unique_centers);

        // Get the index of points closest to the unique_centers
        // Get the index of points closest to the unique_centers
        let medoids: Vec<usize> = unique_centers
            .iter()
            .map(|&center| {
                data.iter()
                    .enumerate()
                    .min_by(|(_, &a), (_, &b)| {
                        (a - center).abs().partial_cmp(&(b - center).abs()).unwrap()
                    })
                    .map(|(idx, _)| idx)
                    .unwrap() // Safe because data is not empty (checked at start)
            })
            .collect();

        MeanShiftResult {
            cluster_centers: unique_centers,
            labels,
            medoids,
        }
    }

    /// Generate seeds using binning approach
    fn get_bin_seeds(&self, data: &[f64], bandwidth: f64, min_bin_freq: usize) -> Vec<f64> {
        if data.is_empty() {
            return vec![];
        }

        // Find data range
        let min_val = data.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_val = data.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        // Create bins
        let bin_size = bandwidth;
        let num_bins = ((max_val - min_val) / bin_size).ceil() as usize + 1;
        if num_bins == 0 {
            return vec![];
        }
        let mut bins: Vec<Vec<f64>> = vec![Vec::new(); num_bins];

        // Assign points to bins
        for &point in data {
            let bin_idx = ((point - min_val) / bin_size) as usize;
            let bin_idx = bin_idx.min(num_bins - 1);
            bins[bin_idx].push(point);
        }

        // Create seeds from bin centers that have enough points
        let mut seeds = Vec::new();
        for (i, bin) in bins.iter().enumerate() {
            if bin.len() >= min_bin_freq {
                let bin_center = min_val + (i as f64 + 0.5) * bin_size;
                seeds.push(bin_center);
            }
        }

        // If no seeds found, use original approach
        if seeds.is_empty() {
            let mut all_seeds: Vec<f64> = data.to_vec();
            all_seeds.sort_by(|a, b| a.partial_cmp(b).unwrap());
            all_seeds.dedup_by(|a, b| (*a - *b).abs() < self.tol);
            return all_seeds;
        }

        seeds
    }

    /// Apply mean shift algorithm to a single seed point
    fn mean_shift_single_seed(&self, data: &[f64], seed: f64, bandwidth: f64) -> f64 {
        let mut current = seed;

        for _ in 0..self.max_iter {
            let mut numerator = 0.0;
            let mut denominator = 0.0;

            // Calculate weighted mean within bandwidth
            for &point in data {
                let distance = (current - point).abs();
                if distance <= bandwidth {
                    // Using flat kernel (uniform kernel)
                    let weight = 1.0;
                    numerator += weight * point;
                    denominator += weight;
                }
            }

            if denominator == 0.0 {
                break;
            }

            let new_center = numerator / denominator;

            // Check for convergence
            if (new_center - current).abs() < self.tol {
                break;
            }

            current = new_center;
        }

        current
    }

    /// Assign cluster labels to data points
    fn assign_labels(&self, data: &[f64], centers: &[f64]) -> Vec<usize> {
        data.iter()
            .map(|&point| {
                let mut best_center = 0;
                let mut min_distance = f64::INFINITY;

                for (i, &center) in centers.iter().enumerate() {
                    let distance = (point - center).abs();
                    if distance < min_distance {
                        min_distance = distance;
                        best_center = i;
                    }
                }

                best_center
            })
            .collect()
    }
}

/// Result of Mean Shift clustering
#[derive(Debug)]
pub struct MeanShiftResult {
    pub cluster_centers: Vec<f64>,
    pub labels: Vec<usize>,
    pub medoids: Vec<usize>,
}

/// Estimate bandwidth for Mean Shift algorithm
///
/// This function estimates the bandwidth by computing pairwise distances
/// and returning the specified quantile of these distances.
/// This version more closely matches scikit-learn's implementation.
///
/// # Arguments
/// * `data` - Input 1D data points
/// * `quantile` - Quantile to use for bandwidth estimation (0.0 to 1.0)
/// * `n_samples` - Number of samples to use for estimation (None = use all)
/// * `random_state` - Random seed for sampling
pub fn estimate_bandwidth(
    data: &[f64],
    quantile: f64,
    n_samples: Option<usize>,
    random_state: u64,
) -> f64 {
    if data.len() < 2 {
        return 1.0; // Default bandwidth for single point
    }

    // Use subset of data if n_samples is specified
    let sample_data: Vec<f64> = match n_samples {
        Some(n) if n < data.len() => {
            // Random sampling with seed (simplified version)
            use std::collections::HashSet;
            let mut rng_state = random_state;
            let mut indices = HashSet::new();

            while indices.len() < n {
                // Simple LCG for reproducible randomness
                rng_state = rng_state.wrapping_mul(1103515245).wrapping_add(12345);
                let idx = (rng_state as usize) % data.len();
                indices.insert(idx);
            }

            indices.into_iter().map(|i| data[i]).collect()
        }
        _ => data.to_vec(),
    };

    // Compute all pairwise distances (this is the expensive O(n²) part)
    let mut distances = Vec::new();

    for i in 0..sample_data.len() {
        for j in (i + 1)..sample_data.len() {
            distances.push((sample_data[i] - sample_data[j]).abs());
        }
    }

    if distances.is_empty() {
        return 1.0;
    }

    // Sort distances and find quantile
    distances.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // More precise quantile calculation
    let index = (distances.len() as f64 * quantile) as usize;
    let index = if index >= distances.len() {
        distances.len() - 1
    } else {
        index
    };

    distances[index]
}

impl Default for MeanShift {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_estimate_bandwidth() {
        let data = vec![1.0, 2.0, 3.0, 10.0, 11.0, 12.0];
        let bandwidth = estimate_bandwidth(&data, 0.5, None, 0);
        assert!(bandwidth > 0.0);
        println!("Estimated bandwidth: {}", bandwidth);
    }

    #[test]
    fn test_meanshift_simple() {
        let data = vec![1.0, 1.1, 1.2, 5.0, 5.1, 5.2, 9.0, 9.1, 9.2];
        let mut ms = MeanShift::new();
        let result = ms.fit(&data);

        println!("Cluster centers: {:?}", result.cluster_centers);
        println!("Labels: {:?}", result.labels);

        // Should find approximately 3 clusters
        assert!(result.cluster_centers.len() >= 2);
        assert_eq!(result.labels.len(), data.len());
    }

    #[test]
    fn test_meanshift_with_bandwidth() {
        let data = vec![1.0, 2.0, 3.0, 10.0, 11.0, 12.0];
        let mut ms = MeanShift::with_bandwidth(2.0);
        let result = ms.fit(&data);

        println!("With bandwidth 2.0 - Centers: {:?}", result.cluster_centers);
        println!("With bandwidth 2.0 - Labels: {:?}", result.labels);

        assert!(!result.cluster_centers.is_empty());
        assert_eq!(result.labels.len(), data.len());
    }

    #[test]
    fn test_single_cluster() {
        let data = vec![1.0, 1.01, 0.99, 1.02, 0.98];
        let mut ms = MeanShift::with_bandwidth(1.0);
        let result = ms.fit(&data);

        println!("Single cluster - Centers: {:?}", result.cluster_centers);
        println!("Single cluster - Labels: {:?}", result.labels);

        // Should find one cluster
        assert_eq!(result.cluster_centers.len(), 1);
        assert!(result.labels.iter().all(|&label| label == 0));
    }

    #[test]
    fn test_mean_shift() {
        let data = [
            5.0, 8.0, 19.0, 6.0, 4.0, 12.0, 9.0, 4.0, 21.0, 8.0, 8.0, 4.0, 10.0, 3.0, 19.0, 10.0,
            20.0, 19.0, 2.0, 20.0,
        ];
        let ans = [0, 1, 2, 0, 0, 1, 1, 0, 2, 1, 1, 0, 1, 0, 2, 1, 2, 2, 0, 2];

        let mut ms = MeanShift::new();
        let result = ms.fit(&data);
        println!("Triple cluster - Centers: {:?}", result.cluster_centers);
        println!("Triple cluster - Labels: {:?}", result.labels);
        println!("Triple cluster - ANSlab: {:?}", ans);
        assert_eq!(result.labels, ans);
        assert_eq!(result.cluster_centers.len(), 3);

        let bandwidth = estimate_bandwidth(&data, 0.3, None, 0);
        println!("Triple Bandwidth : {:?}", bandwidth);
        assert_eq!(bandwidth, 3.0);
    }
}

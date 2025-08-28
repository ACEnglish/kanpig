use crate::kplib::polycluster::PolyCluParams;

/// Mean Shift clustering algorithm for 1D data
pub struct MeanShift {
    pub bandwidth: Option<f64>,
    pub max_iter: usize,
    /// tolerance for convergence
    pub tol: f64,
    /// minimum number of reads per-cluster
    pub min_bin_freq: usize,
}

impl MeanShift {
    /// Create a new MeanShift instance
    pub fn new(params: &PolyCluParams) -> Self {
        Self {
            bandwidth: params.bandwidth,
            max_iter: 300,
            tol: 1e-3,
            min_bin_freq: params.msmin,
        }
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
        Self::new(&PolyCluParams::default())
    }
}

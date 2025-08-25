use ndarray::Array2;

/// Estimate bandwidth like sklearn (for 1D data only).
fn estimate_bandwidth(data: &Array2<f64>, quantile: f64) -> f64 {
    let n = data.nrows();
    let mut dists = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let d = (data[[i, 0]] - data[[j, 0]]).abs();
            dists.push(d);
        }
    }
    dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let idx = (quantile * dists.len() as f64).round() as usize;
    dists[idx.min(dists.len() - 1)]
}

/// Lightweight MeanShift for 1D data.
fn mean_shift_fit(
    data: &Array2<f64>,
    bandwidth: f64,
    min_bin_freq: usize,
) -> (Vec<usize>, usize) {
    let mut points: Vec<f64> = data.column(0).to_vec();
    let mut centers = Vec::new();

    // iterate each point to find its convergence point
    for &mut x in points.iter_mut() {
        let mut prev = f64::NAN;
        let mut cur = x;
        // shift until convergence
        while (cur - prev).abs() > 1e-3 {
            prev = cur;
            // collect neighbors within bandwidth
            let mut neighbors = Vec::new();
            for &p in data.column(0).iter() {
                if (p - cur).abs() <= bandwidth {
                    neighbors.push(p);
                }
            }
            if neighbors.is_empty() {
                break;
            }
            cur = neighbors.iter().sum::<f64>() / neighbors.len() as f64;
        }
        centers.push(cur);
    }

    // Merge close centers (binning)
    let mut unique_centers: Vec<f64> = Vec::new();
    let mut labels = vec![0; centers.len()];

    for (i, &c) in centers.iter().enumerate() {
        let mut found = false;
        for (label, &uc) in unique_centers.iter().enumerate() {
            if (c - uc).abs() < bandwidth / 2.0 {
                labels[i] = label;
                found = true;
                break;
            }
        }
        if !found {
            labels[i] = unique_centers.len();
            unique_centers.push(c);
        }
    }

    // remove clusters below min_bin_freq
    let mut counts = vec![0; unique_centers.len()];
    for &l in &labels {
        counts[l] += 1;
    }
    let mut valid = vec![true; unique_centers.len()];
    for (i, &count) in counts.iter().enumerate() {
        if count < min_bin_freq {
            valid[i] = false;
        }
    }

    // reassign labels (drop invalid clusters → -1)
    for l in labels.iter_mut() {
        if !valid[*l] {
            *l = usize::MAX; // mark as noise
        }
    }

    let n_clusters = valid.iter().filter(|&&x| x).count();

    (labels, n_clusters)
}

/// MeanShift on haplotypes' sizes
/// Returns labels (equal to haps.len) and the n_clusters
pub fn meanshift(
    haps: &[Haplotype],
    bandwidth: Option<f64>,
    min_bin_freq: usize,
    quantile: f64,
) -> (Vec<usize>, usize) {
    let sizes: Vec<f64> = haps.iter().map(|h| h.size as f64).collect();
    let x = Array2::from_shape_vec((sizes.len(), 1), sizes).unwrap();

    let bw = match bandwidth {
        Some(b) => b,
        None => estimate_bandwidth(&x, quantile),
    };

    let (labels, n_clusters) = mean_shift_fit(&x, bw, min_bin_freq);

    (labels, n_clusters)
}


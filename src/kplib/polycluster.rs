use ndarray::{Array, Array2, Axis};
use rand::SeedableRng;

use crate::kplib::{cluster::ClusterResult, hp_sorter, meanshift::MeanShift, Haplotype, PathScore};

// Put this trait on TrioCommand and Mosaic Command so we contain the copying
pub trait ToPolyCluParams {
    fn to_polyclu_params(&self) -> PolyCluParams;
}

// Parameters to pass around the polyclust methods
pub struct PolyCluParams {
    /// Minimum number of reads in a cluster
    pub msmin: usize,

    /// Max clusters
    pub maxclust: usize,

    /// Clustering weight for haplotagged reads (off=0.0, full=1.0)
    pub hps_weight: f32,

    /// Clustering weight for haplotype lengths (off=0.0, full=1.0)
    pub len_weight: f32,

    /// Only cluster on haplotype lengths
    pub lengthonly: bool,

    /// BP difference between MeanShift clusters
    pub bandwidth: Option<f64>,
}

impl Default for PolyCluParams {
    fn default() -> Self {
        PolyCluParams {
            msmin: 5,
            maxclust: 5,
            hps_weight: 0.25,
            len_weight: 0.25,
            lengthonly: false,
            bandwidth: None,
        }
    }
}

pub fn cluster_quality(k: usize, quality_scores: &[f64], labels: &[usize]) -> Vec<f64> {
    // Find the maximum label to determine the size needed
    let mut scores = vec![0.0; k + 1]; // +2 because we need indices 0 through max_label+1

    // Reference cluster is always 0 and quality of 1
    scores[0] = 1.0;

    // Get unique labels and calculate mean for each cluster
    let mut unique_labels: Vec<usize> = labels.to_vec();
    unique_labels.sort_unstable();
    unique_labels.dedup();

    for cluster in unique_labels {
        let cluster_scores: Vec<f64> = quality_scores
            .iter()
            .zip(labels.iter())
            .filter(|(_, &label)| label == cluster)
            .map(|(&score, _)| score)
            .collect();

        if !cluster_scores.is_empty() {
            let mean = cluster_scores.iter().sum::<f64>() / cluster_scores.len() as f64;
            scores[cluster + 1] = mean;
        }
    }

    scores
}

/// Count reads in each cluster
pub fn count_reads(
    k: usize,
    ref_coverage: &[usize],
    assignments: &[usize],
    haplos: &[Haplotype],
) -> Array2<usize> {
    let mut read_counts = Array::<usize, _>::zeros((k + 1, ref_coverage.len()));
    for (i, &coverage) in ref_coverage.iter().enumerate() {
        read_counts[[0, i]] = coverage;
    }
    for (&label, hap) in assignments.iter().zip(haplos.iter()) {
        let cluster_idx = label + 1;
        let sample_idx = hap.meta.samples_flag.trailing_zeros() as usize;
        read_counts[[cluster_idx, sample_idx]] += 1;
    }
    read_counts
}

pub fn top_n_rows_by_sum(arr: &Array2<usize>, n: usize) -> Vec<usize> {
    // Calculate row sums and pair with indices
    let mut row_sums_with_indices: Vec<(usize, usize)> = arr
        .axis_iter(Axis(0)) // Iterate over rows
        .enumerate()
        .skip(1)
        .map(|(idx, rows)| (idx - 1, rows.sum()))
        .collect();

    // Sort by sum in descending order
    row_sums_with_indices.sort_by(|a, b| b.1.cmp(&a.1));

    // Take the top N indices
    row_sums_with_indices
        .into_iter()
        .take(n)
        .map(|(idx, _sum)| idx)
        .collect()
}

// Perform MeanShift and K-medoid clustering
pub fn perform_clustering(
    haplos: &[Haplotype],
    m_args: &PolyCluParams,
    n_samps: usize,
) -> ClusterResult {
    // MeanShift to determine K
    if haplos.is_empty() {
        return ClusterResult {
            assignments: vec![],
            quality: vec![],
            k: 0,
            medoids: vec![],
        };
    }
    if haplos.len() == 1 {
        return ClusterResult {
            assignments: vec![0],
            quality: vec![0.0],
            k: 1,
            medoids: vec![0],
        };
    }

    let sizes: Vec<f64> = haplos.iter().map(|x| x.size as f64).collect();
    let mut ms = MeanShift::new(m_args);
    let ms_result = ms.fit(&sizes);

    // TODO: Experimental: try to make at most 2 like the regular GT does
    let k = ms_result.cluster_centers.len();
    let (mut medoids, k) = if k == 1 {
        // Single center, we can't trust the medoids?
        let medoids = kmedoids::random_initialization(
            haplos.len(),
            2, //m_args.maxclust, // K
            &mut rand::rngs::StdRng::seed_from_u64(21),
        );
        (medoids, 2) //m_args.maxclust)
    } else if k > m_args.maxclust {
        // Only collect the highest covered medoids if MSk > maxclust
        let read_counts = count_reads(k, &vec![0; n_samps], &ms_result.labels, haplos);

        let top = top_n_rows_by_sum(&read_counts, m_args.maxclust);
        let medoids: Vec<usize> = top.iter().map(|&i| ms_result.medoids[i]).collect();
        // Filter haplos to only those in the remaining medoids -- maybe not
        /*let haplos: Vec<Haplotype> = haplos
        .iter()
        .zip(ms_result.labels.iter())
        .filter(|&(_hap, lab)| medoids.contains(lab))
        .map(|(hap, _lab)| hap.clone())
        .collect();*/
        (medoids, m_args.maxclust)
    } else {
        (ms_result.medoids.clone(), k)
    };

    let (assignments, quality) = if m_args.lengthonly {
        // TODO: this is broken. doesn't respect maxclust
        let m_k = ms_result.labels.len();
        (ms_result.labels, vec![1.0; m_k + 1])
    } else {
        // Kmedoid Clustering
        let dist: Array2<f32> = Array2::from_shape_fn((haplos.len(), haplos.len()), |(i, j)| {
            let mut dist: f32 = 1.0 - haplos[i].kmers.fine_similarity(&haplos[j].kmers);
            // if same sample and different hp, hps_weight penalty
            let i_samp = haplos[i].meta.samples_flag;
            let j_samp = haplos[j].meta.samples_flag;
            if i_samp == j_samp {
                match (
                    haplos[i].meta.hp[i_samp.trailing_zeros() as usize],
                    haplos[j].meta.hp[j_samp.trailing_zeros() as usize],
                ) {
                    (Some(group_i), Some(group_j)) if group_i != group_j => {
                        dist *= 1.0 + m_args.hps_weight;
                    }
                    _ => (),
                }
            }
            if ms_result.labels[i] != ms_result.labels[j] {
                dist *= 1.0 + m_args.len_weight;
            }
            dist
        });

        let (_loss, assignments, _, _): (f32, _, _, _) =
            kmedoids::fasterpam(&dist.view(), &mut medoids, 100);
        let (_, quality): (f64, Vec<f64>) = kmedoids::medoid_silhouette(&dist, &medoids, true);
        let quality = cluster_quality(k, &quality, &assignments);
        (assignments, quality)
    };

    // Need to ensure the sizeimilarity of each assignment to the medoid is okay.
    // And if not, what do you assign it to? I don't have a 'drop this' option.
    // So I could make a ClusterResult: No, because assignments are the index.
    // I guess I can remake assignments to be Option<
    ClusterResult {
        assignments,
        quality,
        k,
        medoids,
    }
}

pub fn separate_paths_by_sample(paths: Vec<PathScore>, n_samples: usize) -> Vec<Vec<PathScore>> {
    let mut separated_paths: Vec<Vec<PathScore>> = vec![Vec::new(); n_samples];

    for path in paths {
        for (bit, s_paths) in separated_paths.iter_mut().enumerate().take(n_samples) {
            if (path.meta.samples_flag & (1 << bit)) != 0 {
                let mut p = path.clone();
                p.meta.samples_flag = bit;
                s_paths.push(p);
            }
        }
    }

    // Sort by HP
    separated_paths.iter_mut().for_each(|bin| {
        bin.sort_by(|a, b| {
            hp_sorter(
                &a.meta.hp[a.meta.samples_flag],
                &b.meta.hp[b.meta.samples_flag],
            )
        });
    });
    separated_paths
}

use ndarray::{Array, Array2, Axis};

use crate::kplib::{hp_sorter, metrics, Haplotype, MeanShift, PathScore};

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

    /// Minimum K Freq for seq_to_kmer
    pub minkfreq: u64,
}

impl Default for PolyCluParams {
    fn default() -> Self {
        PolyCluParams {
            msmin: 5,
            maxclust: 5,
            hps_weight: 0.25,
            len_weight: 0.25,
            lengthonly: false,
            minkfreq: 1,
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
        .map(|(idx, row)| (idx - 1, row.sum()))
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

// Structure to hold clustering results
pub struct ClusterResult {
    pub assignments: Vec<usize>,
    pub quality: Vec<f64>,
    pub k: usize,
    pub medoids: Vec<usize>,
}

// Perform MeanShift and K-medoid clustering
pub fn perform_clustering(
    haplos: &[Haplotype],
    m_args: &PolyCluParams,
    n_samps: usize,
) -> ClusterResult {
    // MeanShift to determine K
    let sizes: Vec<f64> = haplos.iter().map(|x| x.size as f64).collect();
    let mut ms = MeanShift::new().min_size(m_args.msmin);
    let ms_result = ms.fit(&sizes);

    let k = ms_result.cluster_centers.len();
    let (mut medoids, k) = if k > m_args.maxclust {
        // Only collect the highest covered medoids if MSk > maxclust
        let read_counts = count_reads(k, &vec![0; n_samps], &ms_result.labels, haplos);

        let top = top_n_rows_by_sum(&read_counts, m_args.maxclust);
        let new_meds = ms_result.medoids.clone();
        (top.iter().map(|&i| new_meds[i]).collect(), m_args.maxclust)
    } else {
        (ms_result.medoids.clone(), k)
    };
    debug!("Setting K to {:?}", k);

    let (assignments, quality) = if m_args.lengthonly {
        // TODO: this is broken. doesn't respect maxclust
        let m_k = ms_result.labels.len();
        (ms_result.labels, vec![1.0; m_k + 1])
    } else {
        // Kmedoid Clustering
        let dist: Array2<f32> = Array2::from_shape_fn((haplos.len(), haplos.len()), |(i, j)| {
            let mut dist: f32 =
                1.0 - (metrics::seqsim(&haplos[i].kfeat, &haplos[j].kfeat, m_args.minkfreq as f32));
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

    ClusterResult {
        assignments,
        quality,
        k,
        medoids,
    }
}

// Process clustered haplotypes and assign reads
pub fn process_clustered_haplotypes(
    cluster_result: ClusterResult,
    haplos: Vec<Haplotype>,
    gts: [[usize; 2]; 3],
) -> Vec<Haplotype> {
    let mut clustered_haps: Vec<Haplotype> = cluster_result
        .medoids
        .iter()
        .map(|i| haplos[*i].clear_clone())
        .collect();

    let mut hp_cnt = Array::<u16, _>::zeros((cluster_result.k, 3, 2));

    // Collapse haplotypes into clusters
    cluster_result
        .assignments
        .into_iter()
        .zip(haplos)
        .for_each(|(cluster_idx, m_hap)| {
            // Sample index inside the HaplotypeMeta
            let idx = m_hap.meta.samples_flag.trailing_zeros() as usize;
            // Only apply reads to the clustered_hap if it goes together
            if gts[idx].contains(&(cluster_idx + 1)) {
                let k_hap = &mut clustered_haps[cluster_idx];
                k_hap.meta.coverage[idx] += 1;
                k_hap.meta.ps[idx] = k_hap.meta.ps[idx].or(m_hap.meta.ps[idx]);
                k_hap.meta.hp[idx] = k_hap.meta.hp[idx].or(m_hap.meta.hp[idx]);
                k_hap.meta.samples_flag |= m_hap.meta.samples_flag;

                if let Some(val) = m_hap.meta.hp[idx] {
                    hp_cnt[[cluster_idx, idx, val as usize - 1]] += 1;
                }
            }
            // TODO: Reassignment of reads assigned to unused clusters?
        });

    // Set HP tag to the most common seen in the cluster
    for (i, m_hap) in clustered_haps.iter_mut().enumerate() {
        for j in 0..3 {
            if m_hap.meta.hp[j].is_some() {
                let max_idx: u8 = hp_cnt
                    .slice(ndarray::s![i, j, ..])
                    .iter()
                    .cloned()
                    .enumerate()
                    .max_by_key(|&(_, val)| val)
                    .map(|(idx, _)| idx)
                    .unwrap_or(1)
                    .try_into()
                    .unwrap();
                m_hap.meta.hp[j] = Some(max_idx + 1);
            }
        }
    }

    clustered_haps
}

pub fn separate_paths_by_sample(paths: Vec<PathScore>) -> Vec<Vec<PathScore>> {
    const NUM_SAMPLES: usize = 3;
    let mut separated_paths: Vec<Vec<PathScore>> = vec![Vec::new(); NUM_SAMPLES];

    for path in paths {
        for (bit, s_paths) in separated_paths.iter_mut().enumerate().take(NUM_SAMPLES) {
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

use crate::kplib::{kmeans, metrics, Haplotype, KDParams};
use itertools::{Either, Itertools};
use crate::kplib::kmeans::Cluster;

fn euclidean_distance(a: &Vec<f32>, b: &Vec<f32>) -> f32 {
    a.iter()
        .zip(b.iter())
        .map(|(x1, x2)| (x1 - x2).powi(2))
        .sum::<f32>()
        .sqrt()
}

pub fn calculate_wcss(clusters: &Vec<Cluster>) -> f32 {
    clusters.iter().map(|cluster| {
        cluster.points.iter().map(|point| {
            let distance = euclidean_distance(&cluster.centroid, point);
            distance.powi(2)
        }).sum::<f32>()
    }).sum()
}

fn calinski_harabasz_index(clusters: &[Cluster]) -> f32 {
    let n_clusters = clusters.len();
    let n_points = clusters
        .iter()
        .map(|cluster| cluster.points.len())
        .sum::<usize>();
    let n_features = clusters[0].centroid.len();
    if n_points <= n_clusters {
        return f32::NAN;
    }
    // Calculate centroids of centroids
    let centroid_of_centroids: Vec<f32> = clusters
        .iter()
        .flat_map(|cluster| cluster.centroid.iter())
        .fold(vec![0.0; n_features], |mut acc, &x| {
            acc.iter_mut().for_each(|a| *a += x);
            acc
        })
        .iter()
        .map(|&sum| sum / n_features as f32)
        .collect();

    // Calculate between-cluster dispersion
    let between_cluster_dispersion: f32 =
        clusters
            .iter()
            .filter(|v| !v.points.is_empty())
            .fold(0.0, |acc, cluster| {
                let cluster_size = cluster.points.len() as f32;
                acc + cluster_size * squared_distance(&cluster.centroid, &centroid_of_centroids)
            });

    // Calculate within-cluster dispersion
    let within_cluster_dispersion: f32 =
        clusters
            .iter()
            .filter(|v| !v.points.is_empty())
            .fold(0.0, |acc, cluster| {
                acc + cluster.points.iter().fold(0.0, |acc, point| {
                    acc + squared_distance(point, &cluster.centroid)
                })
            });

    // Within is perfect, so just return between dispersion?
    if within_cluster_dispersion == 0.0 {
        between_cluster_dispersion / (n_clusters - 1) as f32
    } else {
        // Calculate Calinski-Harabasz index
        (between_cluster_dispersion / (n_clusters - 1) as f32)
            / (within_cluster_dispersion / (n_points - n_clusters) as f32)
    }
}

fn squared_distance(point1: &[f32], point2: &[f32]) -> f32 {
    point1
        .iter()
        .zip(point2.iter())
        .map(|(&x, &y)| (x - y).powi(2))
        .sum()
}

/// Simply takes the best-covered haplotype as the representative
pub fn haploid_haplotypes(
    mut m_haps: Vec<Haplotype>,
    coverage: u64,
    params: &KDParams,
) -> Vec<Haplotype> {
    if coverage == 0 || m_haps.is_empty() {
        return vec![Haplotype::blank(params.kmer, coverage)];
    }

    if m_haps.len() == 1 {
        return m_haps;
    }

    m_haps.sort();

    let mut hap = m_haps.pop().unwrap();
    hap.coverage += m_haps.iter().map(|i| i.coverage).sum::<u64>();

    vec![hap]
}

/// Cluster multiple haplotypes together to try and reduce them to at most two haplotypes
/// This is 'actually' the genotyper. Whatever come out of here is mapped to the variants
/// So inaccurate descriptions of the two haplotypes can not produce good genotypes.
pub fn diploid_haplotypes(
    mut m_haps: Vec<Haplotype>,
    coverage: u64,
    params: &KDParams,
) -> Vec<Haplotype> {
    if coverage == 0 || m_haps.is_empty() {
        return vec![];
    };

    // Nothing to cluster
    if m_haps.len() == 1 {
        let hap = m_haps.pop().unwrap();
        return vec![hap];
    }

    let allk = m_haps.iter().map(|i| i.kfeat.clone()).collect::<Vec<_>>();
    let clusts1 = kmeans(&allk, 1);
    let idx1 = calculate_wcss(&clusts1);
    let clusts2 = kmeans(&allk, 2);
    let idx2 = calculate_wcss(&clusts2);
    
    let clusts = if idx1 <= idx2 {
        clusts1
    } else {
        clusts2
    };

    // Separate the two haplotypes and sort
    let (mut hap_a, mut hap_b): (Vec<Haplotype>, Vec<Haplotype>) =
        m_haps.drain(..).enumerate().partition_map(|(idx, hap)| {
            if clusts[0].points_idx.contains(&idx) {
                Either::Left(hap)
            } else {
                Either::Right(hap)
            }
        });

    // Averaging Haplotypes could be done here?
    // However, 'low-quality' long-reads are now ~97% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10070092/)
    // Which is reasonably within the thresholds of size/seqsim's search space
    hap_a.sort();
    hap_b.sort();
    if hap_a.len() < hap_b.len() {
        std::mem::swap(&mut hap_a, &mut hap_b);
    }

    // Guaranteed to have one cluster
    let mut hap2 = hap_a.pop().unwrap();
    hap2.coverage += hap_a.iter().map(|i| i.coverage).sum::<u64>();

    let mut hap1 = if !hap_b.is_empty() {
        let mut hap_t = hap_b.pop().unwrap();
        hap_t.coverage += hap_b.iter().map(|i| i.coverage).sum::<u64>();
        hap_t
    } else {
        // No deduping needed
        return vec![hap2];
    };

    // hap2 should always be the more supported event
    if hap1 > hap2 {
        std::mem::swap(&mut hap1, &mut hap2);
    }

    debug!("Hap1 in {:?}", hap1);
    debug!("Hap2 in {:?}", hap2);

    // First we establish the two possible alt alleles
    // This is a dedup step for when the alt paths are highly similar
    if (hap1.size.signum() == hap2.size.signum())
        && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs()) > params.hapsim
    {
        hap2.coverage += hap1.coverage;
        return vec![hap2];
    };

    // Now we figure out if the we need two alt alleles or not
    // The reason this takes two steps is the above code is just trying to figure out if
    // there's 1 or 2 alts. Now we figure out if its Het/Hom
    let applied_coverage = (hap1.coverage + hap2.coverage) as f64;
    let remaining_coverage = coverage as f64 - applied_coverage;
    match metrics::genotyper(remaining_coverage, applied_coverage) {
        // We need the one higher covered alt
        // and probably should assign this GT as lowq if REF
        metrics::GTstate::Ref | metrics::GTstate::Het => {
            hap2.coverage += hap1.coverage;
            vec![hap2]
        }
        metrics::GTstate::Hom => {
            if hap1.n == 0 {
                // HOMALT
                vec![hap2]
            } else if hap2.n == 0 {
                // Doesn't happen?
                vec![hap1]
            } else {
                // Compound Het
                vec![hap1, hap2]
            }
        }
        _ => panic!("The genotyper can't do this, yet"),
    }
}

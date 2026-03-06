use crate::kplib::{germ_genotyper, metrics, GraphParams, Haplotype};
use ndarray::{Array, Array2};
use rand::SeedableRng;
use std::collections::HashMap;

// Structure to hold clustering results
pub struct ClusterResult {
    pub assignments: Vec<usize>,
    pub quality: Vec<f64>,
    pub k: usize,
    pub medoids: Vec<usize>,
}

/// Given a cluster result, all observed haplotypes, and a vector of observed haplotypes ids
/// per-sample a.k.a. vec of vec
pub fn collapse_haplotypes(
    cluster_result: ClusterResult,
    haplos: Vec<Haplotype>,
    gts: Vec<Vec<usize>>,
) -> Vec<Haplotype> {
    let mut clustered_haps: Vec<Haplotype> = cluster_result
        .medoids
        .iter()
        .enumerate()
        .map(|(idx, i)| haplos[*i].clear_clone(idx + 1))
        .collect();

    let mut hp_cnt = Array::<u16, _>::zeros((cluster_result.k, gts.len(), 2));

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
                // k_hap.combine(m_hap); I'd like to this, but there's some kinda logic
                // Around Only apply reads to the clustered_hap if it goes together I have to
                // consider.. But I can't rember what it is.
                k_hap.meta.coverage[idx] += 1;
                k_hap.meta.ps[idx] = k_hap.meta.ps[idx].or(m_hap.meta.ps[idx]);
                k_hap.meta.hp[idx] = k_hap.meta.hp[idx].or(m_hap.meta.hp[idx]);
                k_hap.meta.samples_flag |= m_hap.meta.samples_flag;
                k_hap.meta.rnames.append(&mut m_hap.meta.rnames.clone());

                if let Some(val) = m_hap.meta.hp[idx] {
                    hp_cnt[[cluster_idx, idx, val as usize - 1]] += 1;
                }
            }
            // TODO: Reassignment of reads assigned to unused clusters?
        });

    // Set HP tag to the most common seen in the cluster
    for (i, m_hap) in clustered_haps.iter_mut().enumerate() {
        for j in 0..gts.len() {
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

    // Determinism
    clustered_haps.sort();
    clustered_haps
}

/// Simply takes the best-covered haplotype as the representative
pub fn haploid_haplotypes(
    mut haps: Vec<Haplotype>,
    coverage: u64,
    sample_idx: usize,
    _params: &GraphParams,
) -> Vec<Haplotype> {
    if coverage == 0 || haps.is_empty() {
        return vec![];
    }

    let cnt = haps.len() as u64;
    if cnt == 1 {
        return haps;
    }

    let mut g_ps = None;
    let mut rnames = Vec::<String>::with_capacity(haps.len());
    let hap_counts: HashMap<Haplotype, usize> =
        haps.drain(..).fold(HashMap::new(), |mut acc, mut hap| {
            g_ps = g_ps.or(hap.meta.ps[sample_idx]);
            rnames.append(&mut hap.meta.rnames);
            *acc.entry(hap).or_insert(0) += 1;
            acc
        });

    let (mut most_common_hap, _) = hap_counts
        .into_iter()
        .max_by(|(hap1, count1), (hap2, count2)| count1.cmp(count2).then_with(|| hap1.cmp(hap2)))
        .expect("Must be >1 hap to get here");

    most_common_hap.meta.coverage[sample_idx] = cnt;
    most_common_hap.meta.ps[sample_idx] = g_ps;
    most_common_hap.meta.rnames = rnames;

    vec![most_common_hap]
}

/// Cluster multiple haplotypes together to try and reduce them to at most two haplotypes
/// This is 'actually' the genotyper. Whatever come out of here is mapped to the variants
/// So inaccurate descriptions of the two haplotypes can not produce good genotypes.
pub fn diploid_haplotypes(
    mut haplos: Vec<Haplotype>,
    coverage: u64,
    sample_idx: usize,
    hps_weight: f32,
    hapsim: f32,
    ab: f32,
    params: &GraphParams,
) -> Vec<Haplotype> {
    if coverage == 0 || haplos.len() < params.mincoverage || haplos.len() > params.maxcoverage {
        return vec![];
    };

    // So really, all of this should be inside the germ command
    // Nothing to cluster
    if haplos.len() == 1 {
        let hap = haplos.pop().unwrap();
        return vec![hap];
    }

    let distances: Array2<f32> = Array2::from_shape_fn((haplos.len(), haplos.len()), |(i, j)| {
        // Convert similarity to distance
        let dist = 1.0 - &haplos[i].kmers.fine_similarity(&haplos[j].kmers);
        // Penalize only if both points have defined, different groups
        match (haplos[i].meta.hp[sample_idx], haplos[j].meta.hp[sample_idx]) {
            (Some(group_i), Some(group_j)) if group_i != group_j => dist + hps_weight,
            _ => dist,
        }
    });

    let mut medoids = kmedoids::random_initialization(
        haplos.len(),
        2, // K
        &mut rand::rngs::StdRng::seed_from_u64(21),
    );

    let (loss, assignments, _, _): (f32, _, _, _) =
        kmedoids::fasterpam(&distances.view(), &mut medoids, 100);
    debug!("Loss: {}", loss);

    let results = ClusterResult {
        assignments,
        quality: vec![0.0; haplos.len()],
        k: 2,
        medoids,
    };

    let mut haps = collapse_haplotypes(results, haplos, vec![vec![1, 2]]);

    let mut hap1 = haps.swap_remove(0);
    let mut hap2 = haps.swap_remove(0);

    // Hap2 is always the higher covered allele
    if hap2.meta.coverage[sample_idx] < hap1.meta.coverage[sample_idx] {
        std::mem::swap(&mut hap1, &mut hap2);
    }

    debug!("Hap1 in {:?}", hap1);
    debug!("Hap2 in {:?}", hap2);

    // First we establish the two possible alt alleles
    // This is a dedup step for when the alt paths are highly similar
    if (hap1.size.signum() == hap2.size.signum())
        && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs()) > hapsim
    {
        hap2.meta.combine(&hap1.meta);
        return vec![hap2];
    };

    // Now we figure out if the we need two alt alleles or not
    // The reason this takes two steps is the above code is just trying to figure out if
    // there's 1 or 2 alts. Now we figure out if its Het/Compound Het/Hom
    // We have to remove this, I think. Let the genotyper actually do the genotyping
    let applied_coverage = hap1.meta.coverage[sample_idx] + hap2.meta.coverage[sample_idx];
    let remaining_coverage = coverage - applied_coverage;
    let genotyper = germ_genotyper::Genotyper {
        config: germ_genotyper::GenotyperConfig {
            mode: germ_genotyper::GenotypeMode::Beta,
            ..Default::default()
        },
    };
    let gt = genotyper.genotype(
        remaining_coverage,
        hap1.meta.coverage[sample_idx],
        hap2.meta.coverage[sample_idx],
    );

    match gt.state {
        // We need the one higher covered alt
        germ_genotyper::GTstate::Ref | germ_genotyper::GTstate::Het => {
            hap2.meta.combine(&hap1.meta);

            vec![hap2]
        }
        germ_genotyper::GTstate::Hom => {
            if (hap1.meta.coverage[sample_idx] as f32
                / (remaining_coverage + applied_coverage) as f32)
                < ab
            {
                // the allele balance suggests they're not likely compound het
                // Assume hap1 is just a noisy version of hap2
                hap2.meta.combine(&hap1.meta);
                vec![hap2]
            } else {
                vec![hap1, hap2]
            }
        }
        _ => panic!("The genotyper can't do this, yet"),
    }
}

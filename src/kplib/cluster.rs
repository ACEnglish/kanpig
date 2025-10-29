use crate::kplib::{germ_genotyper, metrics, GraphParams, Haplotype};
use ndarray::Array2;
use rand::SeedableRng;
use std::collections::HashMap;

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
    let hap_counts: HashMap<Haplotype, usize> =
        haps.drain(..).fold(HashMap::new(), |mut acc, hap| {
            g_ps = g_ps.or(hap.meta.ps[sample_idx]);
            *acc.entry(hap).or_insert(0) += 1;
            acc
        });

    let (mut most_common_hap, _) = hap_counts
        .into_iter()
        .max_by(|(hap1, count1), (hap2, count2)| count1.cmp(count2).then_with(|| hap1.cmp(hap2)))
        .expect("Must be >1 hap to get here");
    most_common_hap.meta.coverage[sample_idx] = cnt;
    most_common_hap.meta.ps[sample_idx] = g_ps;

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

    // Nothing to cluster
    if haplos.len() == 1 {
        let hap = haplos.pop().unwrap();
        return vec![hap];
    }

    let distances: Array2<f32> = Array2::from_shape_fn((haplos.len(), haplos.len()), |(i, j)| {
        // Convert similarity to distance
        let dist =
            1.0 - (metrics::seqsim(&haplos[i].kfeat, &haplos[j].kfeat, params.minkfreq as f32));
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

    let mut haps = vec![haplos[medoids[0]].clone(), haplos[medoids[1]].clone()];
    let mut hps_cnt = [[0, 0], [0, 0]];

    assignments
        .into_iter()
        .zip(haplos)
        .for_each(|(idx, m_hap)| {
            let k_hap = &mut haps[idx];
            k_hap.meta.coverage[sample_idx] += 1;
            k_hap.meta.ps[sample_idx] = k_hap.meta.ps[sample_idx].or(m_hap.meta.ps[sample_idx]);

            if let Some(hp) = m_hap.meta.hp[sample_idx] {
                hps_cnt[idx][hp as usize - 1] += 1;
            }
        });

    // HP just takes most common
    for (m_hap, hcnts) in haps.iter_mut().zip(hps_cnt) {
        m_hap.meta.coverage[sample_idx] -= 1; // Correct overcounting above
        m_hap.meta.hp[sample_idx] = if hcnts[0] == 0 && hcnts[1] == 0 {
            None
        } else if hcnts[0] >= hcnts[1] {
            Some(1)
        } else {
            Some(2)
        };
    }

    let mut hap1 = haps.swap_remove(0);
    let mut hap2 = haps.swap_remove(0);

    debug!("Hap1 in {:?}", hap1);
    debug!("Hap2 in {:?}", hap2);

    // Hap2 is always the higher covered allele
    if hap2.meta.coverage[sample_idx] < hap1.meta.coverage[sample_idx] {
        std::mem::swap(&mut hap1, &mut hap2);
    }

    // First we establish the two possible alt alleles
    // This is a dedup step for when the alt paths are highly similar
    if (hap1.size.signum() == hap2.size.signum())
        && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs()) > hapsim
    {
        hap2.meta.coverage[sample_idx] += hap1.meta.coverage[sample_idx];
        return vec![hap2];
    };

    // Need some way to challenge if an outlier is by itself..
    /*if hap1.coverage < 3 {
        hap2.coverage += hap1.coverage;
        return vec![hap2];
    }*/

    // Now we figure out if the we need two alt alleles or not
    // The reason this takes two steps is the above code is just trying to figure out if
    // there's 1 or 2 alts. Now we figure out if its Het/Hom
    let applied_coverage = hap1.meta.coverage[sample_idx] + hap2.meta.coverage[sample_idx];
    let remaining_coverage = coverage - applied_coverage;
    let gt = germ_genotyper::phased_genotyper(
        remaining_coverage,
        hap1.meta.coverage[sample_idx],
        hap2.meta.coverage[sample_idx],
    );
    match gt.state {
        // We need the one higher covered alt
        germ_genotyper::GTstate::Ref | germ_genotyper::GTstate::Het => {
            hap2.meta.coverage[sample_idx] += hap1.meta.coverage[sample_idx];
            vec![hap2]
        }
        germ_genotyper::GTstate::Hom => {
            if (hap1.meta.coverage[sample_idx] as f32
                / (remaining_coverage + applied_coverage) as f32)
                < ab
            {
                // the allele balance suggests they're not likely compound het
                // Assume hap1 is just noise and leave it as reference coverage
                vec![hap2]
            } else {
                vec![hap1, hap2]
            }
        }
        _ => panic!("The genotyper can't do this, yet"),
    }
}

use crate::kplib::{metrics, Haplotype, KDParams};
use ndarray::Array2;
use rand::SeedableRng;
use std::collections::HashMap;

/// Simply takes the best-covered haplotype as the representative
pub fn haploid_haplotypes(
    mut haps: Vec<Haplotype>,
    coverage: u64,
    _params: &KDParams,
) -> Vec<Haplotype> {
    if coverage == 0 || haps.is_empty() {
        return vec![];
    }

    let cnt = haps.len() as u64;
    if cnt == 1 {
        return haps;
    }

    let hap_counts: HashMap<Haplotype, usize> =
        haps.drain(..).fold(HashMap::new(), |mut acc, hap| {
            *acc.entry(hap).or_insert(0) += 1;
            acc
        });

    let (mut most_common_hap, _) = hap_counts
        .into_iter()
        .max_by(|(hap1, count1), (hap2, count2)| count1.cmp(count2).then_with(|| hap1.cmp(hap2)))
        .expect("Must be >1 hap to get here");
    most_common_hap.coverage = cnt;

    vec![most_common_hap]
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

    // Create a distance matrix
    let distance_matrix: Array2<f32> =
        Array2::from_shape_fn((m_haps.len(), m_haps.len()), |(i, j)| {
            // Convert similarity to distance
            let dist =
                1.0 - (metrics::seqsim(&m_haps[i].kfeat, &m_haps[j].kfeat, params.minkfreq as f32));
            // Penalize only if both points have defined, different groups
            match (m_haps[i].hp, m_haps[j].hp) {
                (Some(group_i), Some(group_j)) if group_i != group_j => dist + params.hps_weight,
                _ => dist,
            }
        });

    let mut medoids = kmedoids::random_initialization(
        m_haps.len(),
        2, // K
        &mut rand::rngs::StdRng::seed_from_u64(21),
    );

    let (loss, assignments, _, _): (f32, _, _, _) =
        kmedoids::fasterpam(&distance_matrix.view(), &mut medoids, 100);
    trace!("Loss: {}", loss);

    // grab ps from whoever
    let mut g_ps = None;
    let mut hap1 = m_haps[medoids[0]].clone();
    let mut hap1_hps = HashMap::<u8, u16>::new();

    let mut hap2 = m_haps[medoids[1]].clone();
    let mut hap2_hps = HashMap::<u8, u16>::new();
    for (idx, _hap) in m_haps.iter().enumerate() {
        if idx == medoids[0] || idx == medoids[1] {
            continue;
        }
        match assignments[idx] {
            0 => {
                hap1.coverage += 1;
                if let Some(hp) = m_haps[idx].hp {
                    *hap1_hps.entry(hp).or_default() += 1;
                }
            }
            _ => {
                hap2.coverage += 1; // k=2
                if let Some(hp) = m_haps[idx].hp {
                    *hap2_hps.entry(hp).or_default() += 1;
                }
            }
        }

        if g_ps.is_none() && m_haps[idx].ps.is_some() {
            g_ps = m_haps[idx].ps;
        }
    }
    // Assignment of PS and HP just takes most common
    hap1.ps = g_ps;
    if let Some((key, _value)) = hap1_hps.into_iter().max_by_key(|(_, val)| *val) {
        hap1.hp = Some(key);
    }
    hap2.ps = g_ps;
    if let Some((key, _value)) = hap2_hps.into_iter().max_by_key(|(_, val)| *val) {
        hap2.hp = Some(key);
    }

    trace!("Hap1 in {:?}", hap1);
    trace!("Hap2 in {:?}", hap2);

    // Hap2 is always the higher covered allele
    if hap2.coverage < hap1.coverage {
        std::mem::swap(&mut hap1, &mut hap2);
    }

    // First we establish the two possible alt alleles
    // This is a dedup step for when the alt paths are highly similar
    if (hap1.size.signum() == hap2.size.signum())
        && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs()) > params.hapsim
    {
        hap2.coverage += hap1.coverage;
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
    let applied_coverage = (hap1.coverage + hap2.coverage) as f64;
    let remaining_coverage = coverage as f64 - applied_coverage;
    match metrics::genotyper(remaining_coverage, applied_coverage) {
        // We need the one higher covered alt
        metrics::GTstate::Ref | metrics::GTstate::Het => {
            hap2.coverage += hap1.coverage;
            vec![hap2]
        }
        metrics::GTstate::Hom => vec![hap1, hap2],
        _ => panic!("The genotyper can't do this, yet"),
    }
}

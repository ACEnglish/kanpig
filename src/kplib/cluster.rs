use crate::kplib::{metrics, Haplotype, KDParams};
use ndarray::Array2;

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

    // Create a distance matrix based on your similarity function
    let distance_matrix: Array2<f32> =
        Array2::from_shape_fn((m_haps.len(), m_haps.len()), |(i, j)| {
            // Convert similarity to distance
            1.0 - (metrics::seqsim(&m_haps[i].kfeat, &m_haps[j].kfeat, params.minkfreq as f32))
        });
    let mut medoids = kmedoids::random_initialization(m_haps.len(), 2, &mut rand::thread_rng());

    let (loss, assignments, _, _): (f32, _, _, _) =
        kmedoids::fasterpam(&distance_matrix.view(), &mut medoids, 100);
    trace!("Loss: {}", loss);
    let mut hap1 = m_haps[medoids[0]].clone();
    let mut hap2 = m_haps[medoids[1]].clone();
    for (idx, hap) in m_haps.into_iter().enumerate() {
        if idx != medoids[0] && idx != medoids[1] {
            if assignments[idx] == 0 {
                hap1.coverage += hap.coverage;
            } else {
                hap2.coverage += hap.coverage;
            }
        }
    }

    trace!("Hap1 in {:?}", hap1);
    trace!("Hap2 in {:?}", hap2);

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

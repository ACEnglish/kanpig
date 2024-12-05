use crate::kplib::{kmeans, metrics, Haplotype, KDParams};
use itertools::{Either, Itertools};

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
    debug!("{:?}", m_haps);
    let allk = m_haps.iter().map(|i| i.kfeat.clone()).collect::<Vec<_>>();
    let clusts = kmeans(&allk, 2);

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

    trace!("Hap1 in {:?}", hap1);
    trace!("Hap2 in {:?}", hap2);

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

use crate::kplib::{kmeans, metrics, Haplotype, KDParams};

/// Cluster multiple haplotypes together to try and reduce them to at most two haplotypes
/// This is 'actually' the genotyper. Whatever come out of here is mapped to the variants
/// So inaccurate descriptions of the two haplotypes can not produce good genotypes.
pub fn diploid_haplotypes(
    mut m_haps: Vec<Haplotype>,
    coverage: u64,
    params: &KDParams,
) -> Vec<Haplotype> {
    if coverage == 0 || m_haps.is_empty() {
        return vec![
            Haplotype::blank(params.kmer, coverage),
            Haplotype::blank(params.kmer, coverage),
        ];
    };

    // Nothing to cluster
    if m_haps.len() == 1 {
        let hap2 = m_haps.pop().unwrap();
        let ref_cov = (coverage - hap2.coverage) as f64;
        match metrics::genotyper(ref_cov, hap2.coverage as f64) {
            // And if Ref, should probably be set to lowq
            metrics::GTstate::Ref | metrics::GTstate::Het => {
                return vec![Haplotype::blank(params.kmer, ref_cov as u64), hap2]
            }
            metrics::GTstate::Hom => return vec![hap2.clone(), hap2],
            _ => panic!("The genotyper can't do this, yet"),
        }
    }

    let allk = m_haps.iter().map(|i| i.kfeat.clone()).collect::<Vec<_>>();
    let clusts = kmeans(&allk, 2);

    let mut hap_a = Vec::<Haplotype>::new();
    let mut hap_b = Vec::<Haplotype>::new();

    // Separate the two haplotypes and sort
    for (idx, hap) in m_haps.drain(..).enumerate() {
        if clusts[0].points_idx.contains(&idx) {
            hap_a.push(hap);
        } else {
            hap_b.push(hap);
        }
    }
    // Averaging Haplotypes could be done here.
    // However, 'low-quality' long-reads are now ~97% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10070092/)
    // Which is reasonably within the thresholds of size/seqsim's search space
    hap_a.sort();
    hap_b.sort();

    // Guaranteed to have one cluster
    let mut hap2 = hap_a.pop().unwrap();
    hap2.coverage += hap_a.iter().map(|i| i.coverage).sum::<u64>();

    let mut hap1 = if !hap_b.is_empty() {
        let mut hap_t = hap_b.pop().unwrap();
        hap_t.coverage += hap_b.iter().map(|i| i.coverage).sum::<u64>();
        hap_t
    } else {
        // make a ref haplotype from the remaining coverage
        // No deduping needed
        return vec![
            Haplotype::blank(params.kmer, coverage - hap2.coverage),
            hap2,
        ];
    };

    // hap2 should always be the more supported event
    if hap1 > hap2 {
        std::mem::swap(&mut hap1, &mut hap2);
    }

    debug!("Hap1 in {:?}", hap1);
    debug!("Hap2 in {:?}", hap2);

    // First we establish the two possible alt alleles
    // This is a dedup step for when the alt paths are highly similar
    let (hap1, mut hap2) = if (hap1.size.signum() == hap2.size.signum())
        && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs()) >= params.hapsim
    {
        hap2.coverage += hap1.coverage;
        (Haplotype::blank(params.kmer, 0), hap2)
    } else {
        (hap1, hap2)
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
            vec![
                Haplotype::blank(params.kmer, remaining_coverage as u64),
                hap2,
            ]
        }
        metrics::GTstate::Hom => {
            if hap1.n == 0 {
                // HOMALT
                vec![hap2.clone(), hap2]
            } else {
                // Compound Het
                vec![hap1, hap2]
            }
        }
        _ => panic!("The genotyper can't do this, yet"),
    }
}

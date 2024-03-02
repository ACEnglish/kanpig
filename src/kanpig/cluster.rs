use crate::kanpig::{kmeans, metrics, Haplotype, KDParams};

/// Cluster multiple haplotypes together to try and reduce them to at most two haplotypes
/// This is 'actually' the genotyper. Whatever come out of here is mapped to the variants
/// So inaccurate descriptions of the two haplotypes can not produce good genotypes.
pub fn cluster_haplotypes(
    mut m_haps: Vec<Haplotype>,
    coverage: u64,
    params: &KDParams,
) -> (Haplotype, Haplotype) {
    if coverage == 0 || m_haps.is_empty() {
        return (
            Haplotype::blank(params.kmer, coverage),
            Haplotype::blank(params.kmer, coverage),
        );
    };

    // Nothing to cluster
    if m_haps.len() == 1 {
        let mut hap2 = m_haps.pop().unwrap();
        let ref_cov = (coverage - hap2.coverage) as f64;
        let ret = match metrics::genotyper(ref_cov, hap2.coverage as f64) {
            // And if Ref, should probably be set to lowq
            (metrics::GTstate::Ref | metrics::GTstate::Het, _) => {
                (Haplotype::blank(params.kmer, ref_cov as u64), hap2)
            }
            (metrics::GTstate::Hom, gq) => {
                hap2.gq = gq;
                (hap2.clone(), hap2)
            }
            _ => panic!("The genotyper can't do this, yet"),
        };
        return ret;
    }

    let allk = m_haps.iter().map(|i| i.kfeat.clone()).collect::<Vec<_>>();
    let clusts = kmeans(&allk, 2);

    let mut hap_a = Vec::<Haplotype>::new();
    let mut hap_b = Vec::<Haplotype>::new();

    // Sort them so we can pick the highest covered haplotype
    for (idx, hap) in m_haps.drain(..).enumerate() {
        if clusts[0].points_idx.contains(&idx) {
            // println!("Assigning A {:?}", hap);
            hap_a.push(hap);
        } else {
            // println!("Assigning B {:?}", hap);
            hap_b.push(hap);
        }
    }
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
        Haplotype::blank(params.kmer, coverage - hap2.coverage)
    };

    // We want the second haplotype to be the higher covered if it is a non-REF
    if hap1.n != 0
        && (hap1.coverage > hap2.coverage || (hap1.coverage == hap2.coverage && hap1.n < hap2.n))
    {
        std::mem::swap(&mut hap1, &mut hap2);
    }

    // println!("Working with \n\t{:?}\n\t{:?}", hap1, hap2);
    // println!("Hap coverage: {} - {}:{}", coverage, hap1.coverage, hap2.coverage);

    // First we establish the two possible alt alleles
    // This is a dedup step for when the alt paths are highly similar
    let (mut hap1, mut hap2) = match metrics::genotyper(hap1.coverage as f64, hap2.coverage as f64)
    {
        // if hap1 == ref, return hap1, hap2. else combine hap2 into hap1 and make return Hap::blank, hap2
        (metrics::GTstate::Ref, _) => {
            if hap1.n == 0 {
                // println!("Think hap2 is a better representative");
                hap2.coverage += hap1.coverage;
                (Haplotype::blank(params.kmer, 0), hap2)
            } else {
                // println!("Think hap1 is a better representative");
                hap1.coverage += hap2.coverage;
                (Haplotype::blank(params.kmer, 0), hap1)
            }
        }
        // combine them (into hap2, the higher covered allele) return hap2 twice
        (metrics::GTstate::Hom, _) => {
            if (hap1.size.signum() == hap2.size.signum())
                && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs())
                    > params.sizesim
            {
                // println!("These should be consolidated");
                hap2.coverage += hap1.coverage;
                (Haplotype::blank(params.kmer, 0), hap2)
            } else {
                // println!("Could be a compound het or a fp");
                (hap1, hap2)
            }
        }
        // If they're highly similar, combine and assume it was a 'noisy' HOM. Otherwise compound het
        (metrics::GTstate::Het, _) => {
            if (hap1.size.signum() == hap2.size.signum())
                && metrics::sizesim(hap1.size.unsigned_abs(), hap2.size.unsigned_abs())
                    > params.sizesim
            {
                // println!("Turning into Hom Alt");
                hap2.coverage += hap1.coverage;
                (Haplotype::blank(params.kmer, 0), hap2)
            } else {
                // println!("Compound Het");
                (hap1, hap2)
            }
        }
        _ => panic!("The genotyper can't do this, yet"),
    };
    // Now we figure out if the we need two alt alleles or not
    // The reason this takes two steps is the above code is just trying to figure out if
    // there's 1 or 2 alts. Now we figure out if its Het/Hom
    // So essentially we're checking if coverage compared to hap1.coverage/hap2.coverage puts
    // us in ref/het/hom
    let applied_coverage = (hap1.coverage + hap2.coverage) as f64;
    let remaining_coverage = coverage as f64 - applied_coverage;
    match metrics::genotyper(remaining_coverage, applied_coverage) {
        // We need the one higher covered alt
        // and probably should assign this GT as lowq if REF
        (metrics::GTstate::Ref | metrics::GTstate::Het, gq) => {
            hap2.coverage += hap1.coverage;
            hap2.gq = gq;
            (
                Haplotype::blank(params.kmer, remaining_coverage as u64),
                hap2,
            )
        }
        (metrics::GTstate::Hom, gq) => {
            // remaining coverage (if any) is lost here
            hap2.gq = gq.clone();
            if hap1.n == 0 {
                // HOMALT
                (hap2.clone(), hap2)
            } else {
                // Compound Het
                hap1.gq = gq;
                (hap1, hap2)
            }
        }
        _ => panic!("The genotyper can't do this, yet"),
    }
}

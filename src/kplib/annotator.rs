use crate::kplib::{metrics, PathScore, Ploidy};
use bitflags::bitflags;

use petgraph::graph::NodeIndex;

use noodles_vcf::{self as vcf, record::genotypes::sample::Value};

bitflags! {
    pub struct FiltFlags: u32 {
        const PASS       = 0x0;  // passing
        const GTMISMATCH = 0x1;  // genotype from AD doesn't match path genotype
        const LOWGQ      = 0x2;  // genotype quality below 5
        const LOWCOV     = 0x4;  // coverage below 5
        const LOWSQ      = 0x8;  // sample quality below below 5 (only on non-ref genotypes)
        const LOWALT     = 0x16; // alt coverage below 5 (only on non-ref genotypes)
        const PARTIAL    = 0x32; // best scoring path only used part of the haplotype
    }
}

//Format Integer Type Number = G
type IntG = Vec<Option<i32>>;
pub struct GenotypeAnno {
    pub entry: vcf::Record,
    pub gt: String,
    pub filt: FiltFlags,
    pub sq: i32,
    pub gq: i32,
    pub dp: i32,
    pub ad: IntG,
    pub zs: IntG,
    pub ss: IntG,
    pub gt_state: metrics::GTstate,
}

impl GenotypeAnno {
    pub fn new(
        entry: vcf::Record,
        var_idx: &NodeIndex,
        paths: &[PathScore],
        coverage: u64,
        ploidy: &Ploidy,
    ) -> Self {
        match ploidy {
            Ploidy::Zero => zero(entry, coverage),
            Ploidy::Haploid => haploid(entry, var_idx, paths, coverage),
            _ => diploid(entry, var_idx, paths, coverage),
        }
    }

    // These are tied to VcfWriter.keys
    pub fn make_fields(&self, phase_group: i32) -> Vec<Option<Value>> {
        vec![
            Some(Value::from(self.gt.clone())),
            Some(Value::from(self.filt.bits() as i32)),
            Some(Value::from(self.sq)),
            Some(Value::from(self.gq)),
            Some(Value::from(phase_group)),
            Some(Value::from(self.dp)),
            Some(Value::from(self.ad.clone())),
            Some(Value::from(self.zs.clone())),
            Some(Value::from(self.ss.clone())),
        ]
    }
}

/// For annotating a variant in diploid regions
fn diploid(
    entry: vcf::Record,
    var_idx: &NodeIndex,
    paths: &[PathScore],
    coverage: u64,
) -> GenotypeAnno {
    let path1 = &paths[0];
    let path2 = &paths[1];

    let (gt_str, gt_path, alt_cov) =
        match (path1.path.contains(var_idx), path2.path.contains(var_idx)) {
            (true, true) if path1 != path2 => (
                "1|1",
                metrics::GTstate::Hom,
                (path1.coverage.unwrap() + path2.coverage.unwrap()) as f64,
            ),
            // sometimes I used the same path twice
            (true, true) => ("1|1", metrics::GTstate::Hom, path1.coverage.unwrap() as f64),
            (true, false) => ("1|0", metrics::GTstate::Het, path1.coverage.unwrap() as f64),
            (false, true) => ("0|1", metrics::GTstate::Het, path2.coverage.unwrap() as f64),
            (false, false) if coverage != 0 => ("0|0", metrics::GTstate::Ref, 0.0),
            (false, false) => ("./.", metrics::GTstate::Non, 0.0),
        };

    let ref_cov = (coverage as f64) - alt_cov;
    let gt_obs = metrics::genotyper(ref_cov, alt_cov);
    // we're now assuming that ref/alt are the coverages used for these genotypes. no bueno
    let (gq, sq) = metrics::genotype_quals(ref_cov, alt_cov);

    let ad = vec![Some(ref_cov as i32), Some(alt_cov as i32)];

    let zs = vec![
        Some((path1.sizesim * 100.0) as i32),
        Some((path2.sizesim * 100.0) as i32),
    ];

    let ss = vec![
        Some((path1.seqsim * 100.0) as i32),
        Some((path2.seqsim * 100.0) as i32),
    ];

    let mut filt = FiltFlags::PASS;
    // The genotype from AD doesn't match path genotype
    if gt_obs != gt_path {
        filt |= FiltFlags::GTMISMATCH;
    }
    if gq < 5.0 {
        filt |= FiltFlags::LOWGQ;
    }
    if coverage < 5 {
        filt |= FiltFlags::LOWCOV;
    }
    if gt_path != metrics::GTstate::Ref {
        if sq < 5.0 {
            filt |= FiltFlags::LOWSQ;
        }
        if alt_cov < 5.0 {
            filt |= FiltFlags::LOWALT;
        }
    }
    if (path1.align_pct != 1.0) || (path2.align_pct != 1.0) {
        filt |= FiltFlags::PARTIAL;
    }

    GenotypeAnno {
        entry,
        gt: gt_str.to_string(),
        filt,
        sq: sq.round() as i32,
        gq: gq.round() as i32,
        dp: coverage as i32,
        ad,
        zs,
        ss,
        gt_state: gt_path,
    }
}

/// For annotating a variant in a zero ploidy region
fn zero(entry: vcf::Record, coverage: u64) -> GenotypeAnno {
    GenotypeAnno {
        entry,
        gt: "./.".to_string(),
        filt: FiltFlags::PASS,
        sq: 0,
        gq: 0,
        dp: coverage as i32,
        ad: vec![None],
        zs: vec![None],
        ss: vec![None],
        gt_state: metrics::GTstate::Non,
    }
}

/// For annotating a variant in a one ploidy region
fn haploid(
    entry: vcf::Record,
    var_idx: &NodeIndex,
    paths: &[PathScore],
    coverage: u64,
) -> GenotypeAnno {
    let path1 = &paths[0];
    let (gt_str, gt_path, alt_cov) = match path1.path.contains(var_idx) {
        true => ("1", metrics::GTstate::Hom, path1.coverage.unwrap() as f64),
        // sometimes I used the same path twice
        false if coverage != 0 => ("0", metrics::GTstate::Ref, 0.0),
        false => (".", metrics::GTstate::Non, 0.0),
    };

    let ref_cov = (coverage as f64) - alt_cov;
    // we're now assuming that ref/alt are the coverages used for these genotypes. no bueno
    let (gq, sq) = metrics::genotype_quals(ref_cov, alt_cov);

    let ad = vec![Some(ref_cov as i32), Some(alt_cov as i32)];

    let zs = vec![Some((path1.sizesim * 100.0) as i32)];

    let ss = vec![Some((path1.seqsim * 100.0) as i32)];

    let mut filt = FiltFlags::PASS;
    if gq < 5.0 {
        filt |= FiltFlags::LOWGQ;
    }
    if coverage < 5 {
        filt |= FiltFlags::LOWCOV;
    }
    if gt_path != metrics::GTstate::Ref {
        if sq < 5.0 {
            filt |= FiltFlags::LOWSQ;
        }
        if alt_cov < 5.0 {
            filt |= FiltFlags::LOWALT;
        }
    }
    if path1.align_pct != 1.0 {
        filt |= FiltFlags::PARTIAL;
    }

    GenotypeAnno {
        entry,
        gt: gt_str.to_string(),
        filt,
        sq: sq.round() as i32,
        gq: gq.round() as i32,
        dp: coverage as i32,
        ad,
        zs,
        ss,
        gt_state: gt_path,
    }
}

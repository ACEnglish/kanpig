use crate::kplib::{metrics, PathScore, Ploidy};
use bitflags::bitflags;
use noodles_vcf::{
    variant::record_buf::samples::sample::value::{Array, Value},
    variant::RecordBuf,
};
use petgraph::graph::NodeIndex;

bitflags! {
    /// Flags indicating filtering reasons for genotypes.
    pub struct FiltFlags: u32 {
        const PASS       = 0b00000000;  // Passing
        const GTMISMATCH = 0b00000001;  // Genotype from AD doesn't match path genotype
        const LOWGQ      = 0b00000010;  // Genotype quality below 5
        const LOWCOV     = 0b00000100;  // Coverage below 5
        const LOWSQ      = 0b00001000;  // Sample quality below 5 (non-ref genotypes only)
        const LOWALT     = 0b00010000;  // Alt coverage below 5 (non-ref genotypes only)
        const PARTIAL    = 0b00100000;  // Best scoring path uses only part of the haplotype
    }
}

/// Format integer type number for genotype annotations.
type IntG = Vec<Option<i32>>;

/// Struct representing genotype annotations.
pub struct GenotypeAnno {
    pub entry: RecordBuf,
    pub gt: String,
    pub filt: FiltFlags,
    pub sq: i32,
    pub gq: i32,
    pub ps: Option<u32>,
    pub dp: i32,
    pub ad: IntG,
    pub ks: IntG,
    pub gt_state: metrics::GTstate,
}

impl GenotypeAnno {
    /// Creates a new `GenotypeAnno` instance based on the provided ploidy and parameters.
    pub fn new(
        entry: RecordBuf,
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

    /// Generates fields for the `GenotypeAnno` to match `VcfWriter` keys.
    pub fn make_fields(&self, neigh_group: i32) -> Vec<Option<Value>> {
        vec![
            Some(Value::Genotype(
                self.gt.parse().expect("GT string parsing failed"),
            )),
            Some(Value::Integer(self.filt.bits() as i32)),
            Some(Value::Integer(self.sq)),
            Some(Value::Integer(self.gq)),
            self.ps.map(|ps| Value::Integer(ps as i32)),
            Some(Value::Integer(neigh_group)),
            Some(Value::Integer(self.dp)),
            Some(Value::Array(Array::Integer(self.ad.clone()))),
            Some(Value::Array(Array::Integer(self.ks.clone()))),
        ]
    }
}

/// Helper function for a diploid region annotation.
fn diploid(
    entry: RecordBuf,
    var_idx: &NodeIndex,
    paths: &[PathScore],
    coverage: u64,
) -> GenotypeAnno {
    let handle = match &paths {
        [] => handle_diploid_no_paths(coverage),
        [p] => handle_diploid_single_path(var_idx, p, coverage),
        [p1, p2] => handle_diploid_two_paths(var_idx, p1, p2, coverage),
        _ => panic!("Unexpected number of paths for diploid region"),
    };

    finalize_annotation(entry, handle, paths, coverage)
}

/// Helper for zero ploidy regions.
fn zero(entry: RecordBuf, coverage: u64) -> GenotypeAnno {
    GenotypeAnno {
        entry,
        gt: "./.".to_string(),
        filt: FiltFlags::PASS,
        sq: 0,
        gq: 0,
        ps: None,
        dp: coverage as i32,
        ad: vec![None],
        ks: vec![None],
        gt_state: metrics::GTstate::Non,
    }
}

/// Helper for haploid regions.
fn haploid(
    entry: RecordBuf,
    var_idx: &NodeIndex,
    paths: &[PathScore],
    coverage: u64,
) -> GenotypeAnno {
    if paths.is_empty() {
        let handle = match coverage {
            0 => (".", metrics::GTstate::Non, 0.0, true),
            _ => ("0", metrics::GTstate::Ref, 0.0, true),
        };
        return finalize_annotation(entry, handle, paths, coverage);
    }

    // Case 2: At least one path exists
    let path1 = &paths[0];
    let handle = match path1.path.contains(var_idx) {
        true => (
            "1",
            metrics::GTstate::Hom,
            path1.coverage.unwrap_or(0) as f64,
            true,
        ),
        false if coverage != 0 => ("0", metrics::GTstate::Ref, 0.0, true),
        false => (".", metrics::GTstate::Non, 0.0, true),
    };
    finalize_annotation(entry, handle, paths, coverage)
}

/// GT str, GTstate, alt_cov, is_fulltarget
type HandleReturn<'a> = (&'a str, metrics::GTstate, f64, bool);

fn handle_diploid_no_paths<'a>(coverage: u64) -> HandleReturn<'a> {
    if coverage != 0 {
        ("0|0", metrics::GTstate::Ref, 0.0, true)
    } else {
        ("./.", metrics::GTstate::Non, 0.0, true)
    }
}

fn handle_diploid_single_path<'a>(
    var_idx: &NodeIndex,
    path: &PathScore,
    coverage: u64,
) -> HandleReturn<'a> {
    if !path.path.contains(var_idx) {
        ("0|0", metrics::GTstate::Ref, 0.0, true)
    } else {
        let alt_cov = path.coverage.unwrap() as f64;
        let ref_cov = (coverage as f64) - alt_cov;
        let (genotype, state) = match metrics::genotyper(ref_cov, alt_cov) {
            metrics::GTstate::Ref | metrics::GTstate::Het => {
                let gt = match path.hp {
                    None => "0|1",
                    Some(1) => "0|1",
                    _ => "1|0",
                };
                (gt, metrics::GTstate::Het)
            }
            metrics::GTstate::Hom => ("1|1", metrics::GTstate::Hom),
            _ => panic!("Cannot happen here"),
        };
        (genotype, state, alt_cov, path.full_target)
    }
}

fn handle_diploid_two_paths<'a>(
    var_idx: &NodeIndex,
    path1: &PathScore,
    path2: &PathScore,
    coverage: u64,
) -> HandleReturn<'a> {
    match (path1.path.contains(var_idx), path2.path.contains(var_idx)) {
        (true, true) => (
            "1|1",
            metrics::GTstate::Hom,
            (path1.coverage.unwrap() + path2.coverage.unwrap()) as f64,
            path1.full_target || path2.full_target,
        ),
        (true, false) => (
            "1|0",
            metrics::GTstate::Het,
            path1.coverage.unwrap() as f64,
            path1.full_target,
        ),
        (false, true) => (
            "0|1",
            metrics::GTstate::Het,
            path2.coverage.unwrap() as f64,
            path2.full_target,
        ),
        (false, false) if coverage != 0 => ("0|0", metrics::GTstate::Ref, 0.0, true),
        (false, false) => ("./.", metrics::GTstate::Non, 0.0, true),
    }
}

fn finalize_annotation(
    entry: RecordBuf,
    handle: HandleReturn,
    paths: &[PathScore],
    coverage: u64,
) -> GenotypeAnno {
    let (gt_str, gt_path, alt_cov, full_target) = handle;
    let ref_cov = coverage as f64 - alt_cov;

    let gt_obs = metrics::genotyper(ref_cov, alt_cov);

    // we're now assuming that ref/alt are the coverages used for these genotypes. no bueno
    let (gq, sq) = metrics::genotype_quals(ref_cov, alt_cov);

    let ps = if !paths.is_empty() { paths[0].ps } else { None };

    let ad = vec![Some(ref_cov as i32), Some(alt_cov as i32)];

    let ks: Vec<Option<i32>> = paths
        .iter()
        .map(|p| Some((p.score * 100.0) as i32))
        .collect();

    let mut filt = FiltFlags::PASS;
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

    if !full_target {
        filt |= FiltFlags::PARTIAL;
    }

    GenotypeAnno {
        entry,
        gt: gt_str.to_string(),
        filt,
        sq: sq.round() as i32,
        gq: gq.round() as i32,
        ps,
        dp: coverage as i32,
        ad,
        ks,
        gt_state: gt_path,
    }
}

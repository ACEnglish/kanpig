use crate::kplib::{germ_genotyper, PathScore, Ploidy};
use bitflags::bitflags;
use noodles_vcf::{
    header::record::value::map::format,
    variant::record_buf::samples::sample::value::{Array, Value},
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
        const SOMATIC    = 0b01000000;  // SV was found at non-germline VAF
    }
}

/// Format integer type number for genotype annotations.
type IntG = Vec<Option<i32>>;

/// Struct representing genotype annotations.
pub struct GenotypeAnno {
    pub var_idx: NodeIndex,
    pub gt: String,
    pub filt: FiltFlags,
    pub sq: i32,
    pub gq: i32,
    pub ps: Option<u32>,
    pub dp: i32,
    pub ad: IntG,
    pub ks: IntG,
    pub gt_state: germ_genotyper::GTstate,
}

impl GenotypeAnno {
    /// Creates a new `GenotypeAnno` instance based on the provided ploidy and parameters.
    pub fn new(
        var_idx: &NodeIndex,
        paths: &[PathScore],
        coverage: u64,
        ploidy: &Ploidy,
        neigh_group: u64,
        sample_idx: usize, // For pulling the correct coverage from the PathScore.HaplotypeMeta
    ) -> Self {
        match ploidy {
            Ploidy::Zero => zero(*var_idx, coverage),
            Ploidy::Haploid => haploid(var_idx, paths, coverage, neigh_group, sample_idx),
            _ => diploid(var_idx, paths, coverage, neigh_group, sample_idx),
        }
    }

    /// Generates fields for the `GenotypeAnno` used by `VcfWriter`.
    /// Edits to these must be sync'd with make_fmt_definitions
    pub fn make_fields(&self) -> Vec<Option<Value>> {
        // KS can sometimes be an empty array, so we have to set it to None
        let ks = if self.ks.is_empty() {
            None
        } else {
            Some(Value::Array(Array::Integer(self.ks.clone())))
        };

        vec![
            Some(Value::Genotype(
                self.gt.parse().expect("GT string parsing failed"),
            )),
            Some(Value::Integer(self.filt.bits() as i32)),
            Some(Value::Integer(self.sq)),
            Some(Value::Integer(self.gq)),
            self.ps.map(|ps| Value::Integer(ps as i32)),
            Some(Value::Integer(self.dp)),
            Some(Value::Array(Array::Integer(self.ad.clone()))),
            ks,
        ]
    }

    // Edits to these must be sync'd with make_fields
    #[rustfmt::skip]
    pub fn make_format() -> Vec<(&'static str, format::Number, format::Type, &'static str)> {
        let num1 = format::Number::Count(1);
        vec![
            ("GT", num1, format::Type::String, "Kanpig genotype"),
            ("FT", num1, format::Type::Integer, "Kanpig filter"),
            ("SQ", num1, format::Type::Integer, "Phred quality of being non-ref"),
            ("GQ", num1, format::Type::Integer, "Phred quality of genotype"),
            ("PS", num1, format::Type::Integer, "PhaseSet tag from reads"),
            ("DP", num1, format::Type::Integer, "Coverage over region"),
            ("AD", format::Number::ReferenceAlternateBases, format::Type::Integer, "Ref/Alt coverage"),
            ("KS", format::Number::Unknown, format::Type::Integer, "Kanpig score"), // TODO: Extreme values sometimes
        ]
    }
}

/// Helper function for a diploid region annotation.
fn diploid(
    var_idx: &NodeIndex,
    paths: &[PathScore],
    coverage: u64,
    neigh_group: u64,
    sample_idx: usize,
) -> GenotypeAnno {
    let handle = match &paths {
        [] => handle_diploid_no_paths(coverage),
        [p] => handle_diploid_single_path(var_idx, p, coverage, sample_idx),
        [p1, p2] => handle_diploid_two_paths(var_idx, p1, p2, coverage, sample_idx),
        p => panic!("Unexpected number of paths for diploid region {:?}", p),
    };

    finalize_annotation(handle, paths, coverage, neigh_group, sample_idx, *var_idx)
}

/// Helper for zero ploidy regions.
fn zero(var_idx: NodeIndex, coverage: u64) -> GenotypeAnno {
    GenotypeAnno {
        var_idx,
        gt: "./.".to_string(),
        filt: FiltFlags::PASS,
        sq: 0,
        gq: 0,
        ps: None,
        dp: coverage as i32,
        ad: vec![None],
        ks: vec![None],
        gt_state: germ_genotyper::GTstate::Non,
    }
}

/// Helper for haploid regions.
/// Assumed to have ≤1 Path
fn haploid(
    var_idx: &NodeIndex,
    paths: &[PathScore],
    coverage: u64,
    neigh_group: u64,
    sample_idx: usize,
) -> GenotypeAnno {
    if paths.is_empty() {
        let handle = match coverage {
            0 => (".", germ_genotyper::GTstate::Non, 0, true),
            _ => ("0", germ_genotyper::GTstate::Ref, 0, true),
        };
        return finalize_annotation(handle, paths, coverage, neigh_group, sample_idx, *var_idx);
    }

    let path1 = &paths[0];
    let handle = match path1.path.contains(var_idx) {
        true => (
            "1",
            germ_genotyper::GTstate::Hom,
            path1.meta.coverage[sample_idx],
            true,
        ),
        false if coverage != 0 => ("0", germ_genotyper::GTstate::Ref, 0, true),
        false => (".", germ_genotyper::GTstate::Non, 0, true),
    };
    finalize_annotation(handle, paths, coverage, neigh_group, sample_idx, *var_idx)
}

/// GT str, GTstate, alt_cov, is_fulltarget
type HandleReturn<'a> = (&'a str, germ_genotyper::GTstate, u64, bool);

fn handle_diploid_no_paths<'a>(coverage: u64) -> HandleReturn<'a> {
    if coverage != 0 {
        ("0|0", germ_genotyper::GTstate::Ref, 0, true)
    } else {
        ("./.", germ_genotyper::GTstate::Non, 0, true)
    }
}

fn handle_diploid_single_path<'a>(
    var_idx: &NodeIndex,
    path: &PathScore,
    coverage: u64,
    sample_idx: usize,
) -> HandleReturn<'a> {
    if !path.path.contains(var_idx) {
        ("0|0", germ_genotyper::GTstate::Ref, 0, true)
    } else {
        let alt_cov = path.meta.coverage[sample_idx];
        let ref_cov = coverage - alt_cov;
        let (genotype, state) = match germ_genotyper::genotyper(ref_cov, alt_cov).state {
            germ_genotyper::GTstate::Ref | germ_genotyper::GTstate::Het => {
                let gt = match path.meta.hp[sample_idx] {
                    None => "0|1",
                    Some(1) => "0|1",
                    _ => "1|0",
                };
                (gt, germ_genotyper::GTstate::Het)
            }
            germ_genotyper::GTstate::Hom => ("1|1", germ_genotyper::GTstate::Hom),
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
    sample_idx: usize,
) -> HandleReturn<'a> {
    match (path1.path.contains(var_idx), path2.path.contains(var_idx)) {
        (true, true) => (
            "1|1",
            germ_genotyper::GTstate::Hom,
            (path1.meta.coverage[sample_idx] + path2.meta.coverage[sample_idx]),
            path1.full_target || path2.full_target,
        ),
        (true, false) => (
            "1|0",
            germ_genotyper::GTstate::Het,
            path1.meta.coverage[sample_idx],
            path1.full_target,
        ),
        (false, true) => (
            "0|1",
            germ_genotyper::GTstate::Het,
            path2.meta.coverage[sample_idx],
            path2.full_target,
        ),
        (false, false) if coverage != 0 => ("0|0", germ_genotyper::GTstate::Ref, 0, true),
        (false, false) => ("./.", germ_genotyper::GTstate::Non, 0, true),
    }
}

fn finalize_annotation(
    handle: HandleReturn,
    paths: &[PathScore],
    coverage: u64,
    neigh_group: u64,
    sample_idx: usize,
    var_idx: NodeIndex,
) -> GenotypeAnno {
    let (gt_str, gt_path, alt_cov, full_target) = handle;
    let ref_cov = coverage - alt_cov;

    let gt_obs = germ_genotyper::genotyper(ref_cov, alt_cov);

    // we're now assuming that ref/alt are the coverages used for these genotypes. no bueno
    // Either use haplotagging PS or NE
    let ps = paths
        .first()
        .and_then(|p| p.meta.ps.get(sample_idx).copied())
        .unwrap_or(Some(neigh_group as u32));

    let ad = vec![Some(ref_cov as i32), Some(alt_cov as i32)];

    let ks: Vec<Option<i32>> = paths
        .iter()
        .map(|p| Some((p.score * 100.0) as i32))
        .collect();

    let mut filt = FiltFlags::PASS;
    if gt_obs.state != gt_path {
        filt |= FiltFlags::GTMISMATCH;
    }

    if gt_obs.gq < 5.0 {
        filt |= FiltFlags::LOWGQ;
    }

    if coverage < 5 {
        filt |= FiltFlags::LOWCOV;
    }

    if gt_path != germ_genotyper::GTstate::Ref {
        if gt_obs.sq < 5.0 {
            filt |= FiltFlags::LOWSQ;
        }
        if alt_cov < 5 {
            filt |= FiltFlags::LOWALT;
        }
    }

    if !full_target {
        filt |= FiltFlags::PARTIAL;
    }

    GenotypeAnno {
        var_idx,
        gt: gt_str.to_string(),
        filt,
        sq: gt_obs.sq.round() as i32,
        gq: gt_obs.gq.round() as i32,
        ps,
        dp: coverage as i32,
        ad,
        ks,
        gt_state: gt_path,
    }
}

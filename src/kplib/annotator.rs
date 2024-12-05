use crate::kplib::{metrics, PathScore, Ploidy};
use bitflags::bitflags;

use petgraph::graph::NodeIndex;

use noodles_vcf::{
    variant::record_buf::samples::sample::value::{Array, Value},
    variant::RecordBuf,
};

bitflags! {
    pub struct FiltFlags: u32 {
        const PASS       = 0b00000000;  // passing
        const GTMISMATCH = 0b00000001;  // genotype from AD doesn't match path genotype
        const LOWGQ      = 0b00000010;  // genotype quality below 5
        const LOWCOV     = 0b00000100;  // coverage below 5
        const LOWSQ      = 0b00001000;  // sample quality below below 5 (only on non-ref genotypes)
        const LOWALT     = 0b00010000; // alt coverage below 5 (only on non-ref genotypes)
        const PARTIAL    = 0b00100000; // best scoring path only used part of the haplotype
    }
}

//Format Integer Type Number = G
type IntG = Vec<Option<i32>>;
pub struct GenotypeAnno {
    pub entry: RecordBuf,
    pub gt: String,
    pub filt: FiltFlags,
    pub sq: i32,
    pub gq: i32,
    pub ps: Option<u16>,
    pub dp: i32,
    pub ad: IntG,
    pub ks: IntG,
    pub gt_state: metrics::GTstate,
}

impl GenotypeAnno {
    /// Constructs a new `GenotypeAnno` instance based on the provided parameters.
    ///
    /// # Parameters
    /// - `entry`: A `RecordBuf` representing the variant record.
    /// - `var_idx`: A reference to the node index.
    /// - `paths`: A slice of `PathScore` representing the paths.
    /// - `coverage`: An unsigned 64-bit integer representing the coverage.
    /// - `ploidy`: A reference to the `Ploidy` enum representing the ploidy level.
    ///
    /// # Returns
    /// A `GenotypeAnno` instance initialized based on the provided parameters.
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

    /// Generates fields for the `GenotypeAnno` instance.
    /// These fields correspond to the keys defined in `VcfWriter`.
    ///
    /// # Parameters
    /// - `phase_group`: An integer representing the phase group.
    ///
    /// # Returns
    /// A vector containing optional `Value` instances representing the fields of the `GenotypeAnno`.
    pub fn make_fields(&self, neigh_group: i32) -> Vec<Option<Value>> {
        vec![
            Some(Value::Genotype(
                self.gt.parse().expect("Should have made GT correctly"),
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

/// Tries to figure out the order of the GT allele numbers
//fn get_gt_string(is_het: bool, paths: &[PathScore]) -> str {
// I only need to do this for a single pathscore het
// IF, I sort the paths based on their hp

/// For annotating a variant in diploid regions
fn diploid(
    entry: RecordBuf,
    var_idx: &NodeIndex,
    paths: &[PathScore],
    coverage: u64,
) -> GenotypeAnno {
    let (gt_str, gt_path, alt_cov, full_target) = match paths.len() {
        // hom ref
        0 => {
            if coverage != 0 {
                ("0|0", metrics::GTstate::Ref, 0.0, true)
            } else {
                ("./.", metrics::GTstate::Non, 0.0, true)
            }
        }
        // het or hom
        1 => {
            let m_path = &paths[0];
            if !m_path.path.contains(var_idx) {
                ("0|0", metrics::GTstate::Ref, 0.0, true)
            } else {
                let alt_cov = m_path.coverage.unwrap() as f64;
                let ref_cov = (coverage as f64) - alt_cov;
                let (genotype, state) = match metrics::genotyper(ref_cov, alt_cov) {
                    metrics::GTstate::Ref | metrics::GTstate::Het => {
                        let gt = match m_path.hp {
                            None => "0|1",
                            Some(1) => "0|1",
                            _ => "1|0",
                        };
                        (gt, metrics::GTstate::Het)
                    }
                    metrics::GTstate::Hom => ("1|1", metrics::GTstate::Hom),
                    _ => panic!("Cannot happen here"),
                };
                (genotype, state, alt_cov, m_path.full_target)
            }
        }
        // compound het or the other one
        2 => {
            let path1 = &paths[0];
            let path2 = &paths[1];

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
        _ => panic!("Cannot happen in diploid"),
    };

    let ref_cov = coverage as f64 - alt_cov;
    let gt_obs = metrics::genotyper(ref_cov, alt_cov);

    // we're now assuming that ref/alt are the coverages used for these genotypes. no bueno
    let (gq, sq) = metrics::genotype_quals(ref_cov, alt_cov);

    let ps = paths[0].ps;

    let ad = vec![Some(ref_cov as i32), Some(alt_cov as i32)];

    let ks: Vec<Option<i32>> = paths
        .iter()
        .map(|p| Some((p.score * 100.0) as i32))
        .collect();

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

/// For annotating a variant in a zero ploidy region
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

/// For annotating a variant in a one ploidy region
fn haploid(
    entry: RecordBuf,
    var_idx: &NodeIndex,
    paths: &[PathScore],
    coverage: u64,
) -> GenotypeAnno {
    if paths.is_empty() {
        let (gt_str, gq, gt_state) = match coverage == 0 {
            true => (".", 0, metrics::GTstate::Non),
            false => ("0", 100, metrics::GTstate::Ref),
        };

        let mut filt = FiltFlags::PASS;
        if coverage < 5 {
            filt |= FiltFlags::LOWCOV;
        }

        return GenotypeAnno {
            entry,
            gt: gt_str.to_string(),
            filt,
            sq: 0,
            gq,
            ps: None,
            dp: coverage as i32,
            ad: vec![Some(coverage as i32), Some(0)],
            ks: vec![],
            gt_state,
        };
    }

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

    let ps = paths[0].ps;

    let ad = vec![Some(ref_cov as i32), Some(alt_cov as i32)];

    let ks = vec![Some((path1.score * 100.0) as i32)];

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
    if !path1.full_target {
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

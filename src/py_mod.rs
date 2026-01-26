use pyo3::exceptions::PyIOError;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyType;

use rust_htslib::faidx;

use crate::kplib::{
    germ_genotyper::{GTstate, GenotypeMode, GenotypeResult, Genotyper, GenotyperConfig},
    Haplotype, HaplotypeMeta, ReadParser,
};
use std::{path::PathBuf, str::FromStr};

#[pyfunction]
pub fn cansim(a: &PyAny, b: &PyAny, mink: f32) -> PyResult<f32> {
    let vec_a: Vec<f32> = a.extract()?;
    let vec_b: Vec<f32> = b.extract()?;

    Ok(crate::kplib::metrics::seqsim(&vec_a, &vec_b, mink))
}

/// Wrap the Rust function for Python.
/// Input: `sequence: bytes`, `kmer: int`, `negative: bool`, `maxhom: int`
/// Output: list of floats
#[pyfunction]
fn seq_to_kmer(_py: Python<'_>, sequence: String, kmer: u8, negative: bool) -> PyResult<Vec<f32>> {
    // Convert Python bytes -> Rust &[u8]
    let seq: &[u8] = sequence.as_bytes();

    // Call your existing Rust function
    let result = crate::kplib::seq_to_kmer(seq, kmer, negative);

    Ok(result)
}

#[pyclass(name = "PlupParser", unsendable)]
pub struct PyPlupParser {
    inner: crate::kplib::PlupParser,
}

#[pymethods]
impl PyPlupParser {
    /// Create a new PlupParser from file paths and params
    #[new]
    fn new(
        tbx_path: &str,
        reference_path: &str,
        sample_name: String,
        sample_idx: usize,
        sample_count: usize,
    ) -> PyResult<Self> {
        // Open internals here
        let reference = faidx::Reader::from_path(reference_path)
            .map_err(|e| PyValueError::new_err(format!("Failed to open reference: {}", e)))?;

        let params = crate::kplib::GraphParams::default();

        Ok(PyPlupParser {
            inner: crate::kplib::PlupParser::new(
                tbx_path.into(),
                reference,
                sample_name,
                sample_idx,
                sample_count,
                params,
            ),
        })
    }

    pub fn find_pileups(
        &mut self,
        chrom: &str,
        start: u64,
        end: u64,
    ) -> PyResult<(Vec<PyHaplotype>, u64)> {
        let (haps, coverage) = self.inner.find_pileups(chrom, start, end);

        // Convert to PyHaplotype
        let py_haps = haps.into_iter().map(|h| PyHaplotype { inner: h }).collect();

        Ok((py_haps, coverage))
    }

    /// Example method: expose something simple from PlupParser
    fn sample_name(&self) -> String {
        self.inner.get_sample_name()
    }

    fn sample_idx(&self) -> usize {
        self.inner.get_sample_idx()
    }

    fn sample_count(&self) -> usize {
        self.inner.get_sample_count()
    }
}

#[pyclass(name = "Haplotype", unsendable)]
#[derive(Clone)]
pub struct PyHaplotype {
    inner: Haplotype,
}

#[pymethods]
impl PyHaplotype {
    #[new]
    pub fn new(kfeat: Vec<f32>, size: i64, n: u64, hap_meta: PyHaplotypeMeta) -> Self {
        Self {
            inner: Haplotype::new(kfeat, size, n, hap_meta.inner),
        }
    }

    #[staticmethod]
    pub fn blank(kmer: u8, hap_meta: PyHaplotypeMeta) -> Self {
        Self {
            inner: Haplotype::blank(kmer, hap_meta.inner),
        }
    }

    pub fn add(&mut self, other: &PyHaplotype) {
        self.inner.add(&other.inner);
    }

    pub fn partial_haplotypes(
        &self,
        kmer: u8,
        max_fns: usize,
        max_parts: usize,
    ) -> Vec<PyHaplotype> {
        self.inner
            .partial_haplotypes(kmer, max_fns, max_parts)
            .into_iter()
            .map(|h| PyHaplotype { inner: h })
            .collect()
    }

    // Expose fields as read-only properties
    #[getter]
    pub fn size(&self) -> i64 {
        self.inner.size
    }
    #[getter]
    pub fn n(&self) -> u64 {
        self.inner.n
    }
    #[getter]
    pub fn kfeat(&self) -> Vec<f32> {
        self.inner.kfeat.clone()
    }
    #[getter]
    pub fn parts(&self) -> Vec<(i64, Vec<f32>)> {
        self.inner.parts.clone()
    }
    #[getter]
    pub fn partial(&self) -> usize {
        self.inner.partial
    }
    #[getter]
    pub fn meta(&self) -> PyHaplotypeMeta {
        PyHaplotypeMeta {
            inner: self.inner.meta.clone(),
        }
    }
}

#[pyclass(name = "HaplotypeMeta", unsendable)]
#[derive(Clone)]
pub struct PyHaplotypeMeta {
    pub inner: HaplotypeMeta,
}

#[pymethods]
impl PyHaplotypeMeta {
    #[new]
    pub fn new(sample_idx: usize, num_samples: usize) -> Self {
        Self {
            inner: HaplotypeMeta::new(sample_idx, num_samples),
        }
    }

    pub fn combine(&mut self, other: &PyHaplotypeMeta) {
        self.inner.combine(&other.inner);
    }

    // Expose fields as read-write properties
    #[getter]
    pub fn coverage(&self) -> Vec<u64> {
        self.inner.coverage.clone()
    }
    #[setter]
    pub fn set_coverage(&mut self, val: Vec<u64>) {
        self.inner.coverage = val;
    }

    #[getter]
    pub fn ps(&self) -> Vec<Option<u32>> {
        self.inner.ps.clone()
    }
    #[setter]
    pub fn set_ps(&mut self, val: Vec<Option<u32>>) {
        self.inner.ps = val;
    }

    #[getter]
    pub fn hp(&self) -> Vec<Option<u8>> {
        self.inner.hp.clone()
    }
    #[setter]
    pub fn set_hp(&mut self, val: Vec<Option<u8>>) {
        self.inner.hp = val;
    }

    #[getter]
    pub fn samples_flag(&self) -> usize {
        self.inner.samples_flag
    }
    #[setter]
    pub fn set_samples_flag(&mut self, val: usize) {
        self.inner.samples_flag = val;
    }
}

#[pyclass(name = "GenotypeResult", unsendable)]
pub struct PyGenotypeResult {
    pub inner: GenotypeResult,
}

#[pymethods]
impl PyGenotypeResult {
    #[getter]
    pub fn state(&self) -> String {
        match self.inner.state {
            GTstate::Ref => "REF".to_string(),
            GTstate::Het => "HET".to_string(),
            GTstate::Hom => "HOM".to_string(),
            GTstate::Non => "NON".to_string(),
        }
    }

    #[getter]
    pub fn gq(&self) -> f64 {
        self.inner.gq
    }

    #[getter]
    pub fn sq(&self) -> f64 {
        self.inner.sq
    }
}

#[pyclass(name = "Genotyper", unsendable)]
struct PyGenotyper {
    inner: Genotyper,
}

#[pymethods]
impl PyGenotyper {
    #[new]
    #[args(config_path = "None")]
    fn new(config: Option<PyGenotyperConfig>) -> PyResult<Self> {
        Ok(PyGenotyper {
            inner: Genotyper::from_config(config.unwrap().inner),
        })
    }

    #[staticmethod]
    fn from_config_path(path: String) -> PyResult<Self> {
        let path_buf = PathBuf::from(path);
        let genotyper = Genotyper::from_config_file(Some(path_buf));
        Ok(PyGenotyper { inner: genotyper })
    }

    /// Genotype a variant site based on reference and alternate allele coverage
    ///
    /// Parameters
    /// ----------
    /// ref_cov : int
    ///     Coverage of the reference allele
    /// alt_cov1 : int
    ///     Coverage of the first alternate allele
    /// alt_cov2 : int
    ///     Coverage of the second alternate allele
    ///
    /// Returns
    /// -------
    /// tuple of (str, float, float)
    ///     A tuple containing:
    ///     - genotype: The called genotype as a string ("Ref", "Het", "Hom", or "Non")
    ///     - gq: Genotype quality score
    ///     - sq: Sample quality score
    fn genotype(&self, ref_cov: u64, alt1_cov: u64, alt2_cov: u64) -> PyGenotypeResult {
        PyGenotypeResult {
            inner: self.inner.genotype(ref_cov, alt1_cov, alt2_cov),
        }
    }
}

#[pyclass(name = "GenotyperConfig")]
#[derive(Debug, Clone)]
pub struct PyGenotyperConfig {
    inner: GenotyperConfig,
}

#[pymethods]
impl PyGenotyperConfig {
    /// Python __init__: allow constructing manually from fields
    ///
    /// GenotyperConfig(
    ///     mode: str = "Beta",
    ///     mixture_fractions: Optional[List[float]] = None,
    ///     means: Optional[List[float]] = None,
    ///     precisions: Optional[List[float]] = None,
    ///     calibration_table: Optional[List[Tuple[float,float]]] = None,
    /// )
    #[new]
    fn py_new(
        mode: Option<&str>,
        mixture_fractions: Option<Vec<f64>>,
        means: Option<Vec<f64>>,
        precisions: Option<Vec<f64>>,
        calibration_table: Option<Vec<(f64, f64)>>,
    ) -> PyResult<Self> {
        let mode = mode.unwrap_or("Beta");
        let mode_rs =
            GenotypeMode::from_str(mode).map_err(|_| PyValueError::new_err("Invalid Mode"))?;

        let mut inner = GenotyperConfig::default();
        inner.mode = mode_rs;

        if let Some(v) = mixture_fractions {
            inner.mixture_fractions = v;
        }
        if let Some(v) = means {
            inner.means = v;
        }
        if let Some(v) = precisions {
            inner.precisions = v;
        }
        if let Some(v) = calibration_table {
            inner.calibration_table = v;
        }

        Ok(PyGenotyperConfig { inner })
    }

    /// mode as a string: "Beta" | "Bino" | "Phased"
    #[getter]
    pub fn mode(&self) -> String {
        self.inner.mode.as_str().to_string()
    }

    /// mixture_fractions: List[float]
    #[getter]
    pub fn mixture_fractions(&self) -> Vec<f64> {
        self.inner.mixture_fractions.clone()
    }

    /// means: List[float]
    #[getter]
    pub fn means(&self) -> Vec<f64> {
        self.inner.means.clone()
    }

    /// precisions: List[float]
    #[getter]
    pub fn precisions(&self) -> Vec<f64> {
        self.inner.precisions.clone()
    }

    /// calibration_table: List[Tuple[float, float]]
    #[getter]
    pub fn calibration_table(&self) -> Vec<(f64, f64)> {
        self.inner.calibration_table.clone()
    }

    /// Create from a JSON config file.
    ///
    /// Python: GenotyperConfig.from_config_file(path: str) -> GenotyperConfig
    #[classmethod]
    pub fn from_config_path(_cls: &PyType, path: &str) -> PyResult<Self> {
        let pb = PathBuf::from(path);
        match GenotyperConfig::from_config_file(pb) {
            Ok(cfg) => Ok(PyGenotyperConfig { inner: cfg }),
            Err(e) => Err(PyIOError::new_err(format!(
                "Failed to load config from '{}': {}",
                path, e
            ))),
        }
    }

    /// Create from an optional JSON config file.
    ///
    /// Python: GenotyperConfig.from_optional_config(path: Optional[str]) -> GenotyperConfig
    #[classmethod]
    pub fn from_optional_config(_cls: &PyType, path: Option<&str>) -> PyResult<Self> {
        let opt_pb = path.map(PathBuf::from);
        let cfg = GenotyperConfig::from_optional_config(opt_pb);
        Ok(PyGenotyperConfig { inner: cfg })
    }

    /// Default config.
    ///
    /// Python: GenotyperConfig.default() -> GenotyperConfig
    #[classmethod]
    pub fn default(_cls: &PyType) -> PyResult<Self> {
        Ok(PyGenotyperConfig {
            inner: GenotyperConfig::default(),
        })
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "GenotyperConfig(mode='{}', mixture_fractions={:?}, means={:?}, precisions={:?}, calibration_table={:?})",
            self.inner.mode.as_str(),
            self.inner.mixture_fractions,
            self.inner.means,
            self.inner.precisions,
            self.inner.calibration_table,
        ))
    }
}

/// Define the Python module
#[pymodule]
fn kanpig(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(seq_to_kmer, m)?)?;
    m.add_function(wrap_pyfunction!(cansim, m)?)?;
    // m.add_class::<PyKDParams>()?; Too much overhead to bind
    m.add_class::<PyPlupParser>()?;
    m.add_class::<PyHaplotypeMeta>()?;
    m.add_class::<PyHaplotype>()?;
    m.add_class::<PyGenotyperConfig>()?;
    m.add_class::<PyGenotyper>()?;

    Ok(())
}

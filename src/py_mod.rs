use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

use rust_htslib::faidx;
use rust_htslib::tbx::{self};

use crate::kplib::ReadParser;
use crate::kplib::{Haplotype, HaplotypeMeta};

/// Wrap the Rust function for Python.
/// Input: `sequence: bytes`, `kmer: int`, `negative: bool`, `maxhom: int`
/// Output: list of floats
#[pyfunction]
fn seq_to_kmer(
    py: Python<'_>,
    sequence: &PyBytes,
    kmer: u8,
    negative: bool,
    maxhom: usize,
) -> PyResult<Vec<f32>> {
    // Convert Python bytes -> Rust &[u8]
    let seq: &[u8] = sequence.as_bytes();

    // Call your existing Rust function
    let result = crate::kplib::seq_to_kmer(seq, kmer, negative, maxhom);

    Ok(result)
}

#[pyclass(name = "KDParams", unsendable)]
#[derive(Clone)]
pub struct PyKDParams {
    pub inner: crate::kplib::KDParams,
}

#[pymethods]
impl PyKDParams {
    #[new]
    pub fn new(
        passonly: Option<bool>,
        neighdist: Option<u64>,
        sizemin: Option<u32>,
        sizemax: Option<u32>,
        mapq: Option<u8>,
        mapflag: Option<u16>,
        hps_weight: Option<f32>,
        seqsim: Option<f32>,
        sizesim: Option<f32>,
        hapsim: Option<f32>,
        gpenalty: Option<f32>,
        fpenalty: Option<f32>,
        kmer: Option<u8>,
        minkfreq: Option<u64>,
        maxnodes: Option<usize>,
        maxpaths: Option<u64>,
        pileupmax: Option<usize>,
        fnmax: Option<usize>,
        ab: Option<f32>,
        squish: Option<bool>,
        one_to_one: Option<bool>,
        maxhom: Option<usize>,
    ) -> Self {
        let mut params = crate::kplib::KDParams::default();
        if let Some(v) = passonly {
            params.passonly = v;
        }
        if let Some(v) = neighdist {
            params.neighdist = v;
        }
        if let Some(v) = sizemin {
            params.sizemin = v;
        }
        if let Some(v) = sizemax {
            params.sizemax = v;
        }
        if let Some(v) = mapq {
            params.mapq = v;
        }
        if let Some(v) = mapflag {
            params.mapflag = v;
        }
        if let Some(v) = hps_weight {
            params.hps_weight = v;
        }
        if let Some(v) = seqsim {
            params.seqsim = v;
        }
        if let Some(v) = sizesim {
            params.sizesim = v;
        }
        if let Some(v) = hapsim {
            params.hapsim = v;
        }
        if let Some(v) = gpenalty {
            params.gpenalty = v;
        }
        if let Some(v) = fpenalty {
            params.fpenalty = v;
        }
        if let Some(v) = kmer {
            params.kmer = v;
        }
        if let Some(v) = minkfreq {
            params.minkfreq = v;
        }
        if let Some(v) = maxnodes {
            params.maxnodes = v;
        }
        if let Some(v) = maxpaths {
            params.maxpaths = v;
        }
        if let Some(v) = pileupmax {
            params.pileupmax = v;
        }
        if let Some(v) = fnmax {
            params.fnmax = v;
        }
        if let Some(v) = ab {
            params.ab = v;
        }
        if let Some(v) = squish {
            params.squish = v;
        }
        if let Some(v) = one_to_one {
            params.one_to_one = v;
        }
        if let Some(v) = maxhom {
            params.maxhom = v;
        }

        Self { inner: params }
    }
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

        let params = crate::kplib::KDParams::default();

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

/// Define the Python module
#[pymodule]
fn kanpig(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(seq_to_kmer, m)?)?;
    // m.add_class::<PyKDParams>()?; Too much overhead to bind
    m.add_class::<PyPlupParser>()?;
    m.add_class::<PyHaplotypeMeta>()?;
    m.add_class::<PyHaplotype>()?;

    Ok(())
}

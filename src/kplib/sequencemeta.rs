/// This holds read information that's eventually passed to PathScore
/// And then is used by the GenotypeAnno to fill in FORMAT fields
/// In order to allow reads across samples to talk to one another
/// we need to use Vectors. For a single sample operation, we will be
/// accessing everything simply as attribute[0]. For multi-sample, we'll
/// use e.g. attribute[0] for proband, attribute[1] for mother, etc.
/// Since reads can consolidate into a single haplotype, we use the
/// sample_flag as a shortcut to know what samples contributed to the
/// haplotype. e.g. flag & 1 means this is a proband haplotype
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct SequenceMeta {
    pub id: usize,
    pub coverage: Vec<u64>,
    pub ps: Vec<Option<u32>>,
    pub hp: Vec<Option<u8>>,
    pub samples_flag: usize,
    pub rnames: Vec<String>,
}

impl SequenceMeta {
    pub fn new(sample_idx: usize, num_samples: usize) -> Self {
        let mut coverage = vec![0u64; num_samples];
        coverage[sample_idx] += 1;
        let rnames = vec![];

        SequenceMeta {
            id: 0,
            coverage,
            ps: vec![None; num_samples],
            hp: vec![None; num_samples],
            samples_flag: 2_usize.pow(sample_idx as u32),
            rnames,
        }
    }

    /// New SequenceMeta that doesn't belong to anyone
    pub fn new_blank(num_samples: usize) -> Self {
        let coverage = vec![0u64; num_samples];
        SequenceMeta {
            id: 0,
            coverage,
            ps: vec![None; num_samples],
            hp: vec![None; num_samples],
            samples_flag: 0,
            rnames: vec![],
        }
    }

    /// Sets ONLY for a single sample
    pub fn set_ps(&mut self, ps: Option<u32>) {
        self.ps[self.samples_flag.trailing_zeros() as usize] = ps
    }

    /// Sets ONLY for a single sample
    pub fn set_hp(&mut self, hp: Option<u8>) {
        self.hp[self.samples_flag.trailing_zeros() as usize] = hp
    }

    pub fn combine(&mut self, other: &SequenceMeta) {
        for (self_cov, other_cov) in self.coverage.iter_mut().zip(&other.coverage) {
            *self_cov += other_cov;
        }

        for (self_ps, other_ps) in self.ps.iter_mut().zip(&other.ps) {
            if self_ps.is_none() {
                *self_ps = *other_ps;
            }
        }

        for (self_hp, other_hp) in self.hp.iter_mut().zip(&other.hp) {
            if self_hp.is_none() {
                *self_hp = *other_hp;
            }
        }

        self.rnames.extend(other.rnames.clone());

        self.samples_flag |= other.samples_flag;
    }
}

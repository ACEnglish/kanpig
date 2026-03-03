use std::cmp::Ordering;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::ops::{Add, AddAssign};

#[inline]
fn encode_nuc(nuc: u8) -> u64 {
    match nuc.to_ascii_uppercase() {
        b'A' => 0,
        b'G' => 1,
        b'C' => 2,
        b'T' => 3,
        _ => 0,
    }
}

pub fn seq_to_kmer(sequence: &[u8], kmer: u8, negative: bool) -> Vec<(u64, f32)> {
    let ukmer = kmer as usize;
    let cnt = if negative { -1.0f32 } else { 1.0f32 };

    if sequence.len() < ukmer {
        return Vec::new();
    }

    let mut result: Vec<(u64, f32)> = Vec::with_capacity(sequence.len());

    let mut f_result: u64 = 0;
    for (pos, i) in sequence.iter().take(ukmer).enumerate() {
        let f_nuc = encode_nuc(*i);
        f_result += f_nuc << ((ukmer - pos - 1) * 2);
    }
    result.push((f_result, cnt));

    let mask: u64 = (1 << (2 * (kmer - 1) as usize)) - 1;
    for i in sequence.iter().skip(ukmer) {
        let f_nuc = encode_nuc(*i);
        f_result = ((f_result & mask) << 2) + f_nuc;
        result.push((f_result, cnt));
    }

    result.sort_unstable_by_key(|&(k, _)| k);
    result.dedup_by(|a, b| {
        if a.0 == b.0 {
            b.1 += a.1;
            true
        } else {
            false
        }
    });
    result.retain(|&(_, v)| v != 0.0);
    result
}

pub fn merge_kmers(a: &[(u64, f32)], b: &[(u64, f32)]) -> Vec<(u64, f32)> {
    let mut merged = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].0.cmp(&b[j].0) {
            Ordering::Less => {
                merged.push(a[i]);
                i += 1;
            }
            Ordering::Greater => {
                merged.push(b[j]);
                j += 1;
            }
            Ordering::Equal => {
                merged.push((a[i].0, a[i].1 + b[j].1));
                i += 1;
                j += 1;
            }
        }
    }
    merged.extend_from_slice(&a[i..]);
    merged.extend_from_slice(&b[j..]);
    merged
}

pub fn seqsim(a: &[(u64, f32)], b: &[(u64, f32)]) -> f32 {
    let mut deno: f32 = 0.0;
    let mut neum: f32 = 0.0;

    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        let (x, y) = match a[i].0.cmp(&b[j].0) {
            Ordering::Less => {
                i += 1;
                (a[i - 1].1, 0.0)
            }
            Ordering::Greater => {
                j += 1;
                (0.0, b[j - 1].1)
            }
            Ordering::Equal => {
                i += 1;
                j += 1;
                (a[i - 1].1, b[j - 1].1)
            }
        };
        let total_d = x.abs() + y.abs();
        deno += total_d;
        neum += (x - y).abs();
    }
    for &(_, x) in &a[i..] {
        deno += x.abs();
        neum += x.abs();
    }
    for &(_, y) in &b[j..] {
        deno += y.abs();
        neum += y.abs();
    }

    if deno == 0.0 {
        return 0.0;
    }
    if neum == 0.0 {
        return 1.0;
    }
    1.0 - (neum / deno)
}

#[derive(Debug, Clone)]
pub struct KmerVec {
    pub coarse_kfeat: Vec<(u64, f32)>,
    pub fine_kfeat: Vec<(u64, f32)>,
    pub coarse_k: u8,
    pub fine_k: u8,
}

impl KmerVec {
    pub fn new(sequence: &[u8], coarse_kmer: u8, fine_kmer: u8, negative: bool) -> Self {
        KmerVec {
            coarse_kfeat: seq_to_kmer(sequence, coarse_kmer, negative),
            fine_kfeat: seq_to_kmer(sequence, fine_kmer, negative),
            coarse_k: coarse_kmer,
            fine_k: fine_kmer,
        }
    }

    pub fn blank(coarse_kmer: u8, fine_kmer: u8) -> Self {
        KmerVec {
            coarse_kfeat: Vec::new(),
            fine_kfeat: Vec::new(),
            coarse_k: coarse_kmer,
            fine_k: fine_kmer,
        }
    }

    pub fn same_k(&self, other: &Self) -> bool {
        self.coarse_k == other.coarse_k && self.fine_k == other.fine_k
    }

    pub fn coarse_similarity(&self, other: &Self) -> f32 {
        debug_assert!(
            self.same_k(other),
            "KmerVec k mismatch in coarse_similarity"
        );
        seqsim(&self.coarse_kfeat, &other.coarse_kfeat)
    }

    pub fn fine_similarity(&self, other: &Self) -> f32 {
        debug_assert!(self.same_k(other), "KmerVec k mismatch in fine_similarity");
        seqsim(&self.fine_kfeat, &other.fine_kfeat)
    }
}

// + operator: consumes both, returns new KmerVec
impl Add for KmerVec {
    type Output = KmerVec;
    fn add(self, other: KmerVec) -> KmerVec {
        debug_assert!(self.same_k(&other), "KmerVec k mismatch in add");
        KmerVec {
            coarse_k: self.coarse_k,
            fine_k: self.fine_k,
            coarse_kfeat: merge_kmers(&self.coarse_kfeat, &other.coarse_kfeat),
            fine_kfeat: merge_kmers(&self.fine_kfeat, &other.fine_kfeat),
        }
    }
}

// + operator: works on references, no move
impl Add for &KmerVec {
    type Output = KmerVec;
    fn add(self, other: &KmerVec) -> KmerVec {
        debug_assert!(self.same_k(other), "KmerVec k mismatch in add");
        KmerVec {
            coarse_k: self.coarse_k,
            fine_k: self.fine_k,
            coarse_kfeat: merge_kmers(&self.coarse_kfeat, &other.coarse_kfeat),
            fine_kfeat: merge_kmers(&self.fine_kfeat, &other.fine_kfeat),
        }
    }
}

// += operator: mutates self in place
impl AddAssign<&KmerVec> for KmerVec {
    fn add_assign(&mut self, other: &KmerVec) {
        debug_assert!(self.same_k(other), "KmerVec k mismatch in add_assign");
        self.coarse_kfeat = merge_kmers(&self.coarse_kfeat, &other.coarse_kfeat);
        self.fine_kfeat = merge_kmers(&self.fine_kfeat, &other.fine_kfeat);
    }
}

impl PartialEq for KmerVec {
    fn eq(&self, other: &Self) -> bool {
        self.coarse_k == other.coarse_k
            && self.fine_k == other.fine_k
            && self.coarse_kfeat.len() == other.coarse_kfeat.len()
            && self
                .coarse_kfeat
                .iter()
                .zip(&other.coarse_kfeat)
                .all(|(&(ki, vi), &(kj, vj))| ki == kj && vi as u64 == vj as u64)
    }
}

impl Eq for KmerVec {}

impl Hash for KmerVec {
    fn hash<H: Hasher>(&self, state: &mut H) {
        for &(k, v) in &self.coarse_kfeat {
            k.hash(state);
            v.to_bits().hash(state);
        }
    }
}

impl PartialOrd for KmerVec {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for KmerVec {
    fn cmp(&self, other: &Self) -> Ordering {
        let (mut i, mut j) = (0, 0);
        let (a, b) = (&self.coarse_kfeat, &other.coarse_kfeat);
        while i < a.len() && j < b.len() {
            let (ki, vi, kj, vj) = match a[i].0.cmp(&b[j].0) {
                Ordering::Less => {
                    i += 1;
                    (a[i - 1].0, a[i - 1].1 as u64, a[i - 1].0, 0u64)
                }
                Ordering::Greater => {
                    j += 1;
                    (b[j - 1].0, 0u64, b[j - 1].0, b[j - 1].1 as u64)
                }
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                    (a[i - 1].0, a[i - 1].1 as u64, b[j - 1].0, b[j - 1].1 as u64)
                }
            };
            let ord = ki.cmp(&kj).then(vi.cmp(&vj));
            if ord != Ordering::Equal {
                return ord;
            }
        }
        if i < a.len() {
            return Ordering::Greater;
        }
        if j < b.len() {
            return Ordering::Less;
        }
        Ordering::Equal
    }
}

/// Encodes a nucleotide character into its 2-bit representation.
///
/// # Parameters
/// - `nuc`: A byte representing the nucleotide character to encode.
///
/// # Returns
/// A 64-bit unsigned integer representing the binary encoding of the nucleotide.
/// ```
#[inline]
fn encode_nuc(nuc: u8) -> u32 {
    match nuc.to_ascii_uppercase() {
        b'A' => 0,
        b'G' => 1,
        b'C' => 2,
        b'T' => 3,
        _ => 0,
    }
}

/// Converts a DNA sequence into k-mer counts.
/// The k-mers are counted as either positive or negative counts based on the `negative` flag.
///
/// # Parameters
/// - `sequence`: A slice of bytes representing the DNA sequence.
/// - `kmer`: The length of the k-mers to count.
/// - `negative`: If true, the k-mer counts are negative.
///
/// # Returns
/// A vector of k-mer counts represented as floats.
///
/// # Example
/// ```
/// let sequence = b"ACGTACGTAC";
/// let kmer = 3;
/// let negative = false;
/// let kmer_counts = kanpig::seq_to_kmer(sequence, kmer, negative);
/// assert_eq!(kmer_counts.len(), 64); // Example length for k=3
/// ```
pub fn seq_to_kmer(sequence: &[u8], kmer: u8, negative: bool) -> Vec<(u32, f32)> {
    let ukmer = kmer as usize;
    let cnt = if negative { -1.0 } else { 1.0 };

    // Must be at least one kmer long
    if sequence.len() < ukmer {
        return Vec::new();
    }

    let mut result: Vec<(u32, f32)> = Vec::with_capacity(sequence.len());

    // index of the first kmer
    let mut f_result: u32 = 0;
    for (pos, i) in sequence.iter().take(ukmer).enumerate() {
        let f_nuc = encode_nuc(*i);
        f_result += f_nuc << ((ukmer - pos - 1) * 2);
    }
    result.push((f_result, cnt));

    // rolling sum masks off first nuc and adds the next one
    let mask: u32 = (1 << (2 * (kmer - 1) as u32)) - 1;

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

    result
}

/// Combine two seq_to_kmer results
pub fn merge_kmers(a: &[(u32, f32)], b: &[(u32, f32)]) -> Vec<(u32, f32)> {
    let mut merged = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].0.cmp(&b[j].0) {
            std::cmp::Ordering::Less => {
                merged.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                merged.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn seqsim_dense(a: &[f32], b: &[f32], mink: f32) -> f32 {
        let mut deno: f32 = 0.0;
        let mut neum: f32 = 0.0;
        let mut total_d: f32;

        for (&x, &y) in a.iter().zip(b.iter()) {
            total_d = x.abs() + y.abs();
            if total_d >= mink {
                deno += total_d;
                neum += (x - y).abs();
            }
        }

        if deno == 0.0 {
            return 0.0;
        }

        if neum == 0.0 {
            return 1.0;
        }

        1.0 - (neum / deno)
    }

    fn seq_to_kmer_dense(sequence: &[u8], kmer: u8, negative: bool) -> Vec<f32> {
        let ukmer = kmer as usize;
        let mut kcounts = vec![0f32; 1 << (2 * ukmer)];
        let cnt = if negative { -1.0 } else { 1.0 };

        // Must be at least one kmer long
        if sequence.len() < ukmer {
            return kcounts;
        }

        // index of the first kmer
        let mut f_result: u64 = 0;
        for (pos, i) in sequence.iter().take(ukmer).enumerate() {
            let f_nuc = encode_nuc(*i) as u64;
            f_result += f_nuc << ((ukmer - pos - 1) * 2);
        }

        // We know the vector has a space for every possible f_result
        unsafe {
            *kcounts.get_unchecked_mut(f_result as usize) += cnt;
        }

        // rolling sum masks off first nuc and adds the next one
        let mask: u64 = (1 << (2 * (kmer - 1) as usize)) - 1;

        for i in sequence.iter().skip(ukmer) {
            let f_nuc = encode_nuc(*i) as u64;
            f_result = ((f_result & mask) << 2) + f_nuc;

            unsafe {
                *kcounts.get_unchecked_mut(f_result as usize) += cnt;
            }
        }

        kcounts
    }

    fn sparse_to_dense(sparse: &[(u32, f32)], kmer: u8) -> Vec<f32> {
        let mut dense = vec![0f32; 1 << (2 * kmer as usize)];
        for &(k, v) in sparse {
            dense[k as usize] = v;
        }
        dense
    }

    #[test]
    fn test_kmer_roundtrip() {
        let seq = b"GGACATAGACACATAGATAGACCACAGTAGATTGACACAGTTAGACAGATCCGCCGCCGCCGCC";
        for k in [4u8, 6u8] {
            let dense = seq_to_kmer_dense(seq, k, false);
            let sparse = seq_to_kmer(seq, k, false);
            let roundtripped = sparse_to_dense(&sparse, k);
            assert_eq!(dense, roundtripped, "kmer counts differ at k={}", k);
        }
    }

    #[test]
    fn test_kmer_roundtrip_merge() {
        let seq1 = b"GCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCC";
        let seq2 = b"GCCGCGCCACGCCGCGCCGCGTCGCCCGCGCCGGGCCACGCCGCGCGCGCCGCGCCGCGCCGCC";
        for k in [4u8, 6u8] {
            let mut dense = seq_to_kmer_dense(seq1, k, false);
            let dense2 = seq_to_kmer_dense(seq2, k, false);
            dense
                .iter_mut()
                .zip(dense2.iter())
                .for_each(|(x, y)| *x += y);

            let sparse1 = seq_to_kmer(seq1, k, false);
            let sparse2 = seq_to_kmer(seq2, k, false);
            let sparse = merge_kmers(&sparse1, &sparse2);
            let roundtripped = sparse_to_dense(&sparse, k);
            assert_eq!(dense, roundtripped, "kmer counts differ at k={}", k);
        }
    }

    #[test]
    fn test_seqsim_equivalence() {
        let seq1 = b"GGACATAGACACATAGATAGCCCACAGTAGATTGACACAGTTAGACAGATCCGCCGCCGCCGCC";
        let seq2 = b"GGACATAGACACATATATAGACCACAGTAGATTGACACAGTTAGACAGATCCGCCGCCGCCGCC";
        for k in [4u8, 6u8] {
            let dense1 = seq_to_kmer_dense(seq1, k, true);
            let dense2 = seq_to_kmer_dense(seq2, k, true);
            let sparse1 = seq_to_kmer(seq1, k, true);
            let sparse2 = seq_to_kmer(seq2, k, true);

            let sim_dense = seqsim_dense(&dense1, &dense2, 1.0);
            let sim_sparse = crate::kplib::metrics::seqsim(&sparse1, &sparse2, 1.0);
            assert!(
                (sim_dense - sim_sparse).abs() < 1e-6,
                "seqsim differs at k={}: dense={} sparse={}",
                k,
                sim_dense,
                sim_sparse
            );
        }
    }

    #[test]
    fn test_sparse_dense_seqsim() {
        let seq1 = b"GCCGCCGCCGAT";
        let seq2 = b"GCCGCCGCCGCC";
        let k = 4u8;

        let dense1 = seq_to_kmer_dense(seq1, k, false);
        let dense2 = seq_to_kmer_dense(seq2, k, false);
        let sparse1 = seq_to_kmer(seq1, k, false);
        let sparse2 = seq_to_kmer(seq2, k, false);

        // First verify kmer counts match
        for &(idx, v) in &sparse1 {
            assert_eq!(dense1[idx as usize], v, "kmer {} count mismatch", idx);
        }
        // Verify no nonzero entries were dropped
        for (idx, &v) in dense1.iter().enumerate() {
            if v != 0.0 {
                assert!(
                    sparse1.iter().any(|&(k, _)| k == idx as u32),
                    "kmer {} missing from sparse",
                    idx
                );
            }
        }

        let sim_dense = seqsim_dense(&dense1, &dense2, 1.0);
        let sim_sparse = crate::kplib::metrics::seqsim(&sparse1, &sparse2, 1.0);
        assert!(
            (sim_dense - sim_sparse).abs() < 1e-6,
            "dense={} sparse={}",
            sim_dense,
            sim_sparse
        );
    }
}

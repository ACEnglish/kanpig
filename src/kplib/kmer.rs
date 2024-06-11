/// Encodes a nucleotide character into its 2-bit representation.
///
/// # Parameters
/// - `nuc`: A byte representing the nucleotide character to encode.
///
/// # Returns
/// A 64-bit unsigned integer representing the binary encoding of the nucleotide.
/// ```
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

/// Converts a DNA sequence into k-mer counts.
/// Optionally compresses homopolymers before counting k-mers.
/// The k-mers are counted as either positive or negative counts based on the `negative` flag.
///
/// # Parameters
/// - `sequence`: A slice of bytes representing the DNA sequence.
/// - `kmer`: The length of the k-mers to count.
/// - `negative`: If true, the k-mer counts are negative.
/// - `maxhom`: The maximum length of homopolymers; if non-zero, compresses homopolymers before counting k-mers.
///
/// # Returns
/// A vector of k-mer counts represented as floats.
///
/// # Example
/// ```
/// let sequence = b"ACGTACGTAC";
/// let kmer = 3;
/// let negative = false;
/// let maxhom = 2;
/// let kmer_counts = kanpig::seq_to_kmer(sequence, kmer, negative, maxhom);
/// assert_eq!(kmer_counts.len(), 64); // Example length for k=3
/// ```
pub fn seq_to_kmer(sequence: &[u8], kmer: u8, negative: bool, maxhom: usize) -> Vec<f32> {
    if maxhom != 0 {
        return seq_to_kmer(&compress_homopolymer(sequence, maxhom), kmer, negative, 0);
    }

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
        let f_nuc = encode_nuc(*i);
        f_result += f_nuc << ((ukmer - pos - 1) * 2);
    }

    // We know the vector has a space for every possible f_result
    unsafe {
        *kcounts.get_unchecked_mut(f_result as usize) += cnt;
    }

    // rolling sum masks off first nuc and adds the next one
    let mask: u64 = (1 << (2 * (kmer - 1) as usize)) - 1;

    for i in sequence.iter().skip(ukmer) {
        let f_nuc = encode_nuc(*i);
        f_result = ((f_result & mask) << 2) + f_nuc;

        unsafe {
            *kcounts.get_unchecked_mut(f_result as usize) += cnt;
        }
    }

    kcounts
}

/// Compresses sequences of repeated bytes (homopolymers) in the input vector.
/// Limits the length of any sequence of repeated bytes to `maxspan`.
///
/// # Parameters
/// - `vector`: A slice of bytes to be compressed.
/// - `maxspan`: The maximum allowable length for sequences of repeated bytes.
///
/// # Returns
/// A vector of bytes with homopolymer lengths limited to `maxspan`.
pub fn compress_homopolymer(vector: &[u8], maxspan: usize) -> Vec<u8> {
    let mut result = Vec::new();
    let mut count = 0;
    let mut prev_byte = None;

    for byte in vector {
        match prev_byte {
            Some(prev) if prev == byte => {
                count += 1;
                if count < maxspan {
                    result.push(*byte);
                }
            }
            _ => {
                count = 1;
                result.push(*byte);
            }
        }
        prev_byte = Some(byte);
    }

    result
}

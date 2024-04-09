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

/// Count kmers in a sequence
pub fn seq_to_kmer(sequence: &[u8], kmer: u8, negative: bool) -> Vec<f32> {
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

    for i in sequence[1..(sequence.len() - ukmer + 1)].iter() {
        let f_nuc = encode_nuc(*i);
        f_result = ((f_result & mask) << 2) + f_nuc;

        unsafe {
            *kcounts.get_unchecked_mut(f_result as usize) += cnt;
        }
    }

    kcounts
}

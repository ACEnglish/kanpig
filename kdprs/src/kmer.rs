fn encode_nuc(nuc: u8) -> u64 {
    match nuc {
        b'A' => 0,
        b'G' => 1,
        b'C' => 2,
        b'T' => 3,
        _ => 0,
    }
}

/// Count kmers in a sequence
pub fn seq_to_kmer(sequence: &[u8], kmer: u8) -> Vec<f32> {
    let mut kcounts = vec![0f32; 4_usize.pow(kmer.into())];
    let ukmer = kmer as usize;
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

    kcounts[f_result as usize] += 1.0;

    // rolling sum masks off first nuc and adds the next one
    let mask: u64 = 4u64.pow((kmer - 1).into()) - 1;

    for i in sequence[1..(sequence.len() - ukmer + 1)].iter() {
        let f_nuc = encode_nuc(*i);
        f_result = ((f_result & mask) << 2) + f_nuc;

        kcounts[f_result as usize] += 1.0;
    }

    kcounts
}

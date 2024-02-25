use noodles_vcf::{self as vcf};

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

/// Convert vcf::Record to kfeat
pub fn var_to_kfeat(entry: &vcf::Record, kmer: u8) -> (Vec<f32>, i64) {
    let ref_seq = entry.reference_bases().to_string();
    let alt_seq = entry
        .alternate_bases()
        .first()
        .expect("Can only work on sequence resolved variants")
        .to_string();

    let size = alt_seq.len() as i64 - ref_seq.len() as i64;

    let m_ref = seq_to_kmer(&ref_seq.as_bytes()[1..], kmer);
    let m_alt = seq_to_kmer(&alt_seq.as_bytes()[1..], kmer);

    let m_ret: Vec<_> = m_alt
        .iter()
        .zip(m_ref.iter())
        .map(|(&x, &y)| (x - y))
        .collect();

    (m_ret, size)
}

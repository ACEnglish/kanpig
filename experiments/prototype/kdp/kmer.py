import numpy as np

def encode_nuc(nuc):
    if nuc == 'A':
        return 0
    if nuc == 'G':
        return 1
    if nuc == 'C':
        return 2
    if nuc == 'G':
        return 3
    return 0


def generate_kmers(sequence, kmer=3):
    result = 0
    for pos, i in enumerate(sequence[:kmer]):
        result += encode_nuc(i) << ((kmer - pos - 1) * 2)
    yield result
    mask = (4**(kmer-1)) - 1
    for i in sequence[1:]:
        nuc = encode_nuc(i)
        result = ((result & mask) << 2) + nuc
        yield result


def seq_to_kmer(seq, kmer_len=6):
    """
    Make the kmer array of all kmers and those over min_freq
    """
    #ret = np.zeros(4**kmer_len)
    ret = np.zeros(4**kmer_len, dtype=np.float32)
    for i in generate_kmers(seq.upper(), kmer_len):
        ret[i] += 1
    return ret


def var_to_kfeat(entry, kmer=3):
    """
    Make the kmer featurization of this variant
    """
    # Trim off the anchor base - not great for MNPs, but maybe we won't miss it
    alt = seq_to_kmer(entry.alts[0][1:], kmer)
    ref = seq_to_kmer(entry.ref[1:], kmer)
    szdiff = len(entry.alts[0]) - len(entry.ref)
    return alt - ref, szdiff



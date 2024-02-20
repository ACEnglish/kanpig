import pysam
import truvari
import sys

def one_per_hap(calls):
    """
    There are two calls, they are on separate haplotypes
    """
    a, b = calls
    if a.samples[0]['GT'] == (1,1) or b.samples[0]["GT"] == (1,1):
        return False

    if None in a.samples[0]['GT'] or None in b.samples[0]["GT"]:
        return False

    return a.samples[0]['GT'] != b.samples[0]["GT"]

def similarity(m_bench, calls):
    """
    How similar are the calls
    """
    a, b = calls
    match = m_bench.compare_calls([a], [b])
    return match[0].state

def split_evidence(calls):
    """
    1. All variants have the same genotype
    2. There are multiple variants
    """
    gts = set([_.samples[0]['GT'] for _ in calls])
    if len(gts) != 1:
        return False
    return len(calls) > 1

def make_region(calls, buff):
    """
    Need min/max position of variants
    """
    mstart = sys.maxsize
    mend = 0
    for i in calls:
        mstart = min(mstart, i.start)
        mend = max(mend, i.stop)
    return f"{calls[0].chrom}\t{max(0, mstart - buff)}\t{mend + buff}\t{len(calls)}\n"

if __name__ == '__main__':

    truth_vcf_fn = "../GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz"
    truth_bed_file = "../GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed"

    matcher = truvari.Matcher()
    matcher.params.sizefilt = 20
    matcher.params.sizemin = 20
    matcher.params.sizemax = 50000
    matcher.params.seqsim = 0.95
    matcher.params.chunksize = 500
    
    m_bench = truvari.Bench(matcher)

    base = pysam.VariantFile(truth_vcf_fn)
    region_tree = truvari.build_region_tree(base, includebed=truth_bed_file)

    base_i = truvari.region_filter(base, region_tree)
    chunks = truvari.chunker(matcher, ('base', base_i))

    out_files = {"iso": open("chr20.isolated_truth.bed", 'w'),
                 "het_dis": open("chr20.hets.dissimilar.bed", 'w'),
                 "het_sim": open("chr20.hets.similar.bed", 'w'),
                 "split": open("chr20.split.bed", 'w'),
                 "cpx": open("chr20.cpx.bed", 'w')}

    for c, _ in chunks:
        calls = c['base']
        if calls[0].chrom != 'chr20':
            continue
        out_line = make_region(calls, matcher.params.chunksize)
        if len(calls) == 1:
            out_files['iso'].write(out_line)
        elif len(calls) == 2 and one_per_hap(calls):
            if similarity(m_bench, calls):
                out_files['het_sim'].write(out_line)
            else:
                out_files['het_dis'].write(out_line)
        elif split_evidence(calls):
            out_files['split'].write(out_line)
        else:
            out_files['cpx'].write(out_line)

    for i in out_files.values():
        i.close()



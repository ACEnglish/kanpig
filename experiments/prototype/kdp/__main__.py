import sys
import logging
import itertools

from functools import partial

import pysam
import truvari

import kdp

def unify_chromosomes(tree, bam_fn, ref_fn):
    """
    Remove regions from VCF tree that bam or reference don't hold
    """
    vcf_names = set(tree.keys())
    bam = pysam.AlignmentFile(bam_fn)
    bam_names = set(bam.references)
    ref = pysam.FastaFile(ref_fn)
    ref_names = set(ref.references)
    bam_rm = vcf_names - bam_names
    ref_rm = vcf_names - ref_names
    if bam_rm or ref_rm:
        logging.warning("VCF has %d references not in BAM and %d not in REF. Removing", len(bam_rm), len(ref_rm))
        for i in bam_rm & ref_rm:
            del(tree[i])
    return tree
            
def main():
    io_params, kd_params = kdp.parse_args(sys.argv[1:])
    logging.info("Starting")

    if io_params.vcf is not None:
        base = pysam.VariantFile(io_params.vcf)
    comp = pysam.VariantFile(io_params.input)

    matcher = truvari.Matcher()
    matcher.params.passonly = kd_params.passonly
    matcher.params.sizefilt = kd_params.sizemin
    matcher.params.sizemin = kd_params.sizemin
    matcher.params.sizemax = kd_params.sizemax
    matcher.params.chunksize = kd_params.chunksize

    if io_params.vcf is not None:
        region_tree = truvari.build_region_tree(base, comp, io_params.regions)
    else:
        region_tree = truvari.build_region_tree(comp, includebed=io_params.regions)
        region_tree = unify_chromosomes(region_tree, io_params.bam, io_params.reference)
    truvari.merge_region_tree_overlaps(region_tree)
    comp_i = truvari.region_filter(comp, region_tree)
    if io_params.vcf is not None:
        base_i = truvari.region_filter(base, region_tree)
        chunks = truvari.chunker(matcher, ('base', base_i), ('comp', comp_i))
    else:
        chunks = truvari.chunker(matcher, ('comp', comp_i))
    

    header = comp.header.copy()
    header.add_line(('##FORMAT=<ID=SZ,Number=R,Type=Float,'
                    'Description="Size similarity of phase group">'))
    header.add_line(('##FORMAT=<ID=CS,Number=R,Type=Float,'
                    'Description="Cosine similarity of phase group">'))
    header.add_line(('##FORMAT=<ID=PG,Number=1,Type=String,'
                    'Description="Phase group id from kdp">'))
    header.add_line(('##FORMAT=<ID=AD,Number=R,Type=Integer,'
                     'Description="Allele Depth">'))
    out_name = io_params.output
    if out_name.endswith(".vcf.gz"):
        out_name = io_params.output[:-len(".gz")]
    out_vcf = pysam.VariantFile(out_name, 'w', header=header)

    # switching types.
    io_params.sample = 0 if io_params.sample is None else io_params.sample

    if io_params.vcf is not None:
        task = partial(kdp.kdp_job_vcf,
                       params=kd_params,
                       header=header,
                       sample=io_params.sample)
    else:
        task = partial(kdp.kdp_job_bam,
                       bam=pysam.AlignmentFile(io_params.bam),
                       reference=pysam.FastaFile(io_params.reference),
                       header=header,
                       sample=io_params.sample,
                       params=kd_params)

    for variant in itertools.chain.from_iterable(map(task, chunks)):
        out_vcf.write(variant)
    out_vcf.close()

    if io_params.output.endswith(".gz"):
        truvari.compress_index_vcf(out_name)
    logging.info("Finished")

if __name__ == '__main__':
    main()

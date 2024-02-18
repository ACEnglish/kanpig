import sys
import logging
import argparse
import itertools

from functools import partial

import pysam
import truvari

import kdp

def parse_args(args):
    """
    """
    parser = argparse.ArgumentParser(prog="kdp", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="VCF to genotype")
    parser.add_argument("--vcf", type=str,
                        help="Phased VCF with genotypes to apply")
    parser.add_argument("--bam", type=str,
                        help="Bam file with reads to apply")
    parser.add_argument("-f", "--reference", type=str,
                        help="Reference the bam was aligned to")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output VCF (stdout)")
    parser.add_argument("-r", "--regions", type=str, default=None,
                        help="Bed filename or comma-separated list of chrom:start-end regions to process")
    parser.add_argument("--kmer", type=int, default=3,
                        help="Kmer size(%(default)s)")
    parser.add_argument("-s", "--sample", type=str, default=None,
                        help="Name of sample to apply genotypes (first)")
    parser.add_argument("--passonly", action="store_true",
                        help="Only phase passing variants")
    parser.add_argument("--sizemin", type=int, default=20,
                        help="Minimum variant size (%(default)s)")
    parser.add_argument("--sizemax", type=int, default=50000,
                        help="Maximum variant size (%(default)s)")
    parser.add_argument("--maxpaths", type=int, default=1000,
                        help="Stop region processing after trying N paths (%(default)s)")
    parser.add_argument("--cossim", type=float, default=0.90,
                        help="Minimum cosine similarity (%(default)s)")
    parser.add_argument("--pctsize", type=truvari.restricted_float, default=0.90,
                        help="Minimum size similarity between a path and haplotype (%(default)s)")
    parser.add_argument("--pg", action="store_true",
                        help="Allow multiple phase groups (%(default)s)")
    parser.add_argument("--chunksize", type=truvari.restricted_int, default=100,
                        help="Max reference distance to phase calls (%(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="Verbose logging")

    args = parser.parse_args(args)
    truvari.setup_logging(args.debug)
    if (not args.bam and not args.vcf) or args.bam and args.vcf:
        logging.error("One of --vcf xor --bam must be provided")
        sys.exit(1)
    if args.sizemin < 20:
        logging.warning("--sizemin is recommended to be at least 20")
    if args.bam and not args.reference:
        logging.error("--bam requires --reference")
        sys.exit(1)
    elif args.vcf and args.reference:
        logging.warning("--reference is not used with --vcf (until we start dealing with symbolic alts)")
    return args

def main():
    args = parse_args(sys.argv[1:])
    logging.info("Starting")

    if args.vcf is not None:
        base = pysam.VariantFile(args.vcf)
    comp = pysam.VariantFile(args.input)

    matcher = truvari.Matcher()
    matcher.params.passonly = args.passonly
    matcher.params.sizefilt = args.sizemin
    matcher.params.sizemin = args.sizemin
    matcher.params.sizemax = args.sizemax
    matcher.params.chunksize = args.chunksize

    if args.vcf is not None:
        region_tree = truvari.build_region_tree(base, comp, args.regions)
    else:
        # Probably need to figure out how to check bam contigs
        region_tree = truvari.build_region_tree(comp, includebed=args.regions)
    truvari.merge_region_tree_overlaps(region_tree)
    comp_i = truvari.region_filter(comp, region_tree)
    if args.vcf is not None:
        base_i = truvari.region_filter(base, region_tree)
        chunks = truvari.chunker(matcher, ('base', base_i), ('comp', comp_i))
    else:
        chunks = truvari.chunker(matcher, ('comp', comp_i))
    
    out_name = args.output
    if out_name.endswith(".vcf.gz"):
        out_name = args.output[:-len(".gz")]

    header = comp.header.copy()
    header.add_line(('##FORMAT=<ID=SZ,Number=R,Type=Float,'
                    'Description="Size similarity of phase group">'))
    header.add_line(('##FORMAT=<ID=CS,Number=R,Type=Float,'
                    'Description="Cosine similarity of phase group">'))
    header.add_line(('##FORMAT=<ID=PG,Number=1,Type=String,'
                    'Description="Phase group id from kdp">'))
    header.add_line(('##FORMAT=<ID=AD,Number=R,Type=Integer,'
                     'Description="Allele Depth">'))
    args.sample = 0 if args.sample is None else args.sample
    out_vcf = pysam.VariantFile(out_name, 'w', header=header)

    if args.vcf is not None:
        task = partial(kdp.kdp_job_vcf,
                       max_paths=args.maxpaths,
                       phase_groups=args.pg,
                       header=header,
                       kmer=args.kmer,
                       min_cos=args.cossim,
                       min_size=args.pctsize,
                       sample=args.sample)
    else:
        task = partial(kdp.kdp_job_bam,
                       bam=pysam.AlignmentFile(args.bam),
                       reference=pysam.FastaFile(args.reference),
                       sizemin=args.sizemin,
                       sizemax=args.sizemax,
                       buffer=args.chunksize,
                       max_paths=args.maxpaths,
                       phase_groups=args.pg,
                       header=header,
                       kmer=args.kmer,
                       min_cos=args.cossim,
                       min_size=args.pctsize,
                       sample=args.sample)

    for variant in itertools.chain.from_iterable(map(task, chunks)):
        out_vcf.write(variant)
    out_vcf.close()

    if args.output.endswith(".gz"):
        truvari.compress_index_vcf(out_name)
    logging.info("Finished")

if __name__ == '__main__':
    main()

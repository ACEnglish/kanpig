import argparse
import logging
import truvari
from dataclasses import dataclass

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
    parser.add_argument("-s", "--sample", type=str, default=None,
                        help="Name of sample to apply genotypes (first)")
    parser.add_argument("--kmer", type=int, default=4,
                        help="Kmer size(%(default)s)")
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
    parser.add_argument("--wcoslen", type=truvari.restricted_int, default=2000,
                        help="Size threshold to switch to weighted cossim (%(default)s)")
    parser.add_argument("--n_tries", type=truvari.restricted_int, default=5,
                        help="Number of attempts to find variants (%(default)s)")
    parser.add_argument("--pg", action="store_true",
                        help="Allow multiple phase groups (don't use) (%(default)s)")
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
    if args.kmer >= 8:
        logging.warning("Kmer size above 8 become memory intensive")

    io_params = IOParams(args.input,
                         args.vcf,
                         args.bam,
                         args.reference,
                         args.output,
                         args.regions,
                         args.sample,
                         args.debug)

    kd_params = KDParams(args.kmer,
                         args.passonly,
                         args.sizemin,
                         args.sizemax,
                         args.maxpaths,
                         args.cossim,
                         args.pctsize,
                         args.wcoslen,
                         args.chunksize,
                         args.n_tries)
    return io_params, kd_params

@dataclass
class IOParams():
    """
    Holds the IO parameters for a job
    """
    input: str
    vcf: str = None
    bam: str = None
    reference: str = None
    output: str = None
    regions: str = None
    sample: str = None
    debug: bool = False

@dataclass
class KDParams():
    """
    Hold all the KDP parameters for a job
    """
    kmer: int = 4
    passonly: bool = True
    sizemin: int = 20
    sizemax: int = 50000
    maxpaths: int = 1000
    cossim: float = 0.90
    pctsize: float = 0.90
    wcoslen: int = 500
    chunksize: int = 100
    n_tries: int = 5

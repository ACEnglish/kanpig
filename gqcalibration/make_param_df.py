"""
Extract relevant GT information from a VCF with truth and kanpig genotypes.
First sample must be truth set.
Second sample must be kanpig run on that sample.
"""
import sys
import truvari
import pandas as pd
import numpy as np
import argparse

def parse_args(args):
    parser = argparse.ArgumentParser(prog="make_param_df", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("VCF", type=str,
                        help="Input VCF of genotypes")
    parser.add_argument("OUT", type=str,
                        help="Output CSV file to write")
    parser.add_argument("--bed", default=None, type=str,
                        help="Bed file for subsetting VCF entries to parse")
    return parser.parse_args(args)

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    vcf = truvari.VariantFile(args.VCF)
    m_iter = vcf.fetch_bed(args.bed) if args.bed else vcf
    rows = []
    for entry in m_iter:
        if entry.chrom in ['chrX', 'chrY'] \
            or entry.var_size() > 10000 \
            or entry.is_monrefstar() \
            or None in entry.samples[1]['GT']:
            continue

        b_gt = truvari.get_gt(entry.gt(0))
        o_gt = truvari.get_gt(entry.gt(1))
        rows.append([b_gt == o_gt,
                     b_gt.name,
                     o_gt.name,
                     entry.samples[1]['DP'],
                     *entry.samples[1]['AD'],
                     entry.samples[1]['GQ'],
                     entry.samples[1]['FT'],
                     min(entry.samples[1]['KS']),
                     ])
    out = pd.DataFrame(rows, columns=['state', 'Ogt', 'Mgt', 'DP', 'AD_ref', 'AD_alt', 'GQ', 'FT', 'KS'])
    out.to_csv(args.OUT, index=False)

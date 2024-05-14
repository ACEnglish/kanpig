kanpig - Kmer ANalysis of PIleups for Genotyping
<img src="https://github.com/ACEnglish/kanpig/raw/develop/imgs/icon.png/" style="width:100px;" align="right" style="vertical-align: middle;"> 
------
A fast tool for genotyping structural variants with long-reads.

# Install
```
git clone https://github.com/ACEnglish/kanpig
cd kanpig
cargo build --release
# executable in ./target/release/kanpig
```

# Usage
```
Usage: kanpig [OPTIONS] --input <INPUT> --bam <BAM> --reference <REFERENCE> --out <OUT>

Options:
  -i, --input <INPUT>          VCF to genotype
  -b, --bam <BAM>              Reads to genotype
  -f, --reference <REFERENCE>  Reference bam is aligned to
  -o, --out <OUT>              Output vcf (unsorted)
      --bed <BED>              Regions to analyze
      --sample <SAMPLE>        Sample to apply genotypes to, default first column
      --threads <THREADS>      Number of threads [default: 1]
      --debug                  Verbose logging
      --kmer <KMER>            Kmer size for featurization [default: 4]
      --chunksize <CHUNKSIZE>  Minimum distance between variants to create independent graphs [default: 100]
      --passonly               Only analyze reads with PASS FILTER
      --sizemin <SIZEMIN>      Minimum size of variant to analyze [default: 50]
      --sizemax <SIZEMAX>      Maximum size of variant to analyze [default: 50000]
      --maxpaths <MAXPATHS>    Maximum number of paths in a graph to traverse [default: 1000]
      --seqsim <SEQSIM>        Minimum sequence similarity for paths [default: 0.9]
      --sizesim <SIZESIM>      Minimum size similarity for paths [default: 0.95]
      --minkfreq <MINKFREQ>    Minimum frequency of kmer [default: 1]
      --hapsim <HAPSIM>        Haplotype size similarity collapse threshold [default: 0.95]
      --try-exact              Search for a 1-to-1 match before graph traversal
      --prune                  Prune paths which don't traverse 1-to-1 nodes
      --mapq <MAPQ>            Minimum mapq of reads to consider [default: 5]
      --mapflag <MAPFLAG>      Alignments with flag matching this value are ignored [default: 3840]
      --spanoff                Don't require alignments to span vargraph region
  -h, --help                   Print help
  -V, --version                Print version
```

# Method

For 'chunks' of variants, kanpig will build a directed variant graph. A chunk of variants is comprised of all vcf
entries within `--sizemin` and `--sizemax` that have at most `--chunksize` distance between them. The minimum start
position and maximum end position of a chunk of variants becomes the region of interest. Pileups of all reads which
span the region of interest are then generated. Reads which have an alignment flag which matches `--mapflag` (where
matching is `(alnflag & mapflag) != 0`) or mapq below `--mapq` are ignored. Every read's pileup sequences are then
featurized with a `--kmer` and sum of variant lengths recorded into a 'haplotype'. Haplotypes that exactly match are
then consolidated. If more than one haplotype remains, kmeans clustering is performed to produce up to two clusters.
If there are two haplotypes, their size similarity is calculated and if above `--hapsim` the higher coverage haplotype
is kept and consolidated with the lower coverage haplotype. A genotyper then analyzes the coverage supporting alternate
haplotypes with coverage supporting the reference. This genotyping will produce one pair of hom-ref, het, hom-alt, or
compound-alt haplotypes. These haplotypes are then applied to the variant graph. This is done by traversing all possible
paths through the graph and summing the variants' k-feats. Each path is compared to the haplotype and the highest
scoring path is kept. Scores are the `--seqsim` of the k-feats and `--sizesim` of the allele deltas. All variants inside
the highest scoring paths are then genotyped as being present in that haplotype. Annotations based on the coverage are
also added to the written vcf entry.

# Annotations
* FT - Bit flag for properties of the variant's genotyping. Flags == 0 are considered PASS. The bits definitions are:
  * 0x1 - The genotype observed from variants matching paths is not equal to the genotype observed from measuring the
  proportions of reads supporting the two alleles.
  * 0x2 - The genotype quality is less than 5
  * 0x4 - The DP is less than 5
  * 0x8 - The sample quality is less than 5 (non-Ref)
  * 0x16 - The number of reads supporting the alternate allele less than 5 (non-Ref)
* SQ - Phred scaled likelihood variant alternate is present in the sample
* GQ - Phred scale difference between most and second-most likely genotypes
* PG - Each chunk of variants is assigned a phase group
* DP - Read coverage over the region
* AD - Read coverage supporting the reference and alternate alleles.
* SZ - Size similarity of the two haplotypes to this variant
* SS - Sequence similarity of the two haplotypes to this variant

# Parameter Details

## `--bed`
Sorted bed file restricts kanpig to only analyzing variants with starts and ends within a single bed entry.

## `--ploidy-bed`
This bed file informs kanpig of special regions within chromosomes that should have non-diploid genotypes. For example, a female
human sample shouldn't have any genotypes on chrY. A male human sample should have hemizygous genotypes on chrY and the
non-pseudoautosomal regions of chrX. The `ploidy_beds/` directory has example bed files for GRCh38. All regions not
within the `ploidy-bed` (or if no bed is provided) are assumed to be diploid.

## `--hapsim`
Wh



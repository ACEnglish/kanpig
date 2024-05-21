kanpig - Kmer ANalysis of PIleups for Genotyping
<img src="https://github.com/ACEnglish/kanpig/raw/develop/imgs/icon.png/" style="width:100px;" align="right" style="vertical-align: middle;"> 
------
A fast tool for genotyping structural variants with long-reads.

*Kanpig is currently under active research and development. We make no guarantees about its accuracy or the stability of features 
before version 1.0.*

# Install
```
git clone https://github.com/ACEnglish/kanpig
cd kanpig
cargo build --release
# executable in ./target/release/kanpig
```
Alternatively, binaries are available in releases.

# Usage
```
Usage: kanpig [OPTIONS] --input <INPUT> --bam <BAM> --reference <REFERENCE> --out <OUT>

Options:
  -i, --input <INPUT>            VCF to genotype
  -b, --bam <BAM>                Reads to genotype
  -f, --reference <REFERENCE>    Reference bam is aligned to
  -o, --out <OUT>                Output vcf (unsorted)
      --bed <BED>                Regions to analyze
      --ploidy-bed <PLOIDY_BED>  Bed file of non-diploid regions
      --sample <SAMPLE>          Sample to apply genotypes to, default first column
      --threads <THREADS>        Number of threads [default: 1]
      --debug                    Verbose logging
      --kmer <KMER>              Kmer size for featurization [default: 4]
      --chunksize <CHUNKSIZE>    Minimum distance between variants to create independent graphs [default: 1000]
      --passonly                 Only analyze reads with PASS FILTER
      --sizemin <SIZEMIN>        Minimum size of variant to analyze [default: 50]
      --sizemax <SIZEMAX>        Maximum size of variant to analyze [default: 10000]
      --maxpaths <MAXPATHS>      Maximum number of paths in a graph to traverse [default: 10000]
      --seqsim <SEQSIM>          Minimum sequence similarity for paths [default: 0.9]
      --sizesim <SIZESIM>        Minimum size similarity for paths [default: 0.9]
      --minkfreq <MINKFREQ>      Minimum frequency of kmer [default: 1]
      --hapsim <HAPSIM>          Haplotype size similarity collapse threshold (off=1) [default: 1]
      --try-exact                Search for a 1-to-1 match before graph traversal
      --prune                    Prune paths which don't traverse 1-to-1 nodes
      --mapq <MAPQ>              Minimum mapq of reads to consider [default: 5]
      --mapflag <MAPFLAG>        Alignments with flag matching this value are ignored [default: 3840]
      --spanoff                  Don't require alignments to span vargraph region
      --maxhom <MAXHOM>          Maximum homopolymer length to kmerize (off=0) [default: 0]
  -h, --help                     Print help
  -V, --version                  Print version
```

# Current Limitations
* Kanpig expects sequence resolved SVs. Variants with symbolic alts (e.g. `<DEL>`) and BNDs are not parsed.
* Kanpig only looks at read pileups and does not consider split or soft-clipped alignment information. This means
  variants above ~10kbp should be skipped with the `--sizemax` parameter.
* Please do not publish manuscripts with kanpig results until we've completed our manuscript. We're aiming to have a preprint 
available early Q3 2024.

# Annotations

Kanpig 
* FT - Bit flag for properties of the variant's genotyping. Flags == 0 are considered PASS. The bits definitions are:
  * 0x1 - The genotype observed from variants matching paths is not equal to the genotype observed from measuring the
  proportions of reads supporting the two alleles.
  * 0x2 - The genotype quality is less than 5
  * 0x4 - The depth (DP) is less than 5
  * 0x8 - The sample quality (SQ) is less than 5 (only present on non-ref variants)
  * 0x16 - The number of reads supporting the alternate allele less than 5 (only present on non-ref variants)
  * 0x32 - The best scoring path through the variant graph only used part of the haplotype. This may be indicative of a
    false-negative in the variant graph.
* SQ - Phred scaled likelihood variant alternate is present in the sample
* GQ - Phred scale difference between most and second-most likely genotypes
* PG - Each chunk of variants is assigned a phase group
* DP - Read coverage over the region
* AD - Read coverage supporting the reference and alternate alleles.
* SZ - Size similarity of the two haplotypes to this variant
* SS - Sequence similarity of the two haplotypes to this variant

# Core Parameter Details

The default parameters are tuned to work generally well for genotyping a single sample's VCF, meaning the variants are
all expected to be present in the sample. For a multi-sample VCF (a.k.a. a project-level VCF), the optimal parameters
are still being determined and will likely be dependent on things such as number of samples in the VCF, merging strategy 
of the variants, and sequencing technology.

## `--bed`
Sorted bed file restricts kanpig to only analyzing variants with starts and ends within a single bed entry.

## `--ploidy-bed`
This bed file informs kanpig of special regions within chromosomes that should have non-diploid genotypes. For example, a female
human sample shouldn't have any genotypes on chrY. A male human sample should have hemizygous genotypes on chrY and the
non-pseudoautosomal regions of chrX. The `ploidy_beds/` directory has example bed files for GRCh38. All regions not
within the `ploidy-bed` (or if no bed is provided) are assumed to be diploid.

## `--hapsim`
After performing kmeans clustering on reads to determine the two haplotypes, if the two haplotypes have a size similarity above `hapsim`, they
are consolidated into a homozygous allele.

## `--chunksize`
Kanpig will build local variant graphs from windows of the genome. These windows are determined by `chunksize` where
the maximum end position of an upstream window's variants is at least `chunksize` base-pairs away from the next window's
variants' minimum start position.

This chunksize also determins the region over which read pileups are generated. Only reads which pass the `--mapq` and
`--mapflag` filter are considered. Also, reads must fully span the minimum variant start and maximum variant end. 

`chunksize` is an important parameter because too small of a value may not recruit read pileups which support variants
but are far away. Similarly, too large of a value may create windows with many SVs which are also too large for reads to have a
fully-spanning alignment. 

## `--sizemin` and `--sizemax`
Variant sizes are determined by `INFO/SVLEN`. If `INFO/SVLEN` tag is not in the VCF entry, the variant's size is set as
`abs(length(ALT) - length(REF))`.

## `--sizesim` and `--seqsim`
When applying a haplotype to a variant graph, only paths above these two thresholds are allowed. If there are multiple
paths above the threshold, the one with the higher `(sizesim + seqsim) / 2` is kept. Generally, `0.90` is well balanced
whereas lower thresholds will boost recall at the cost of precision and vice versa for higher thresholds.

## `--maxpaths`
When performing path-finding, this threshold limits the number of paths which are checked. A lower `--maxpaths` will
speed up runtime but may come at a cost of recall. A higher `--maxpaths` is slower and may come at a cost to
specificity.

## `--threads`
Number of analysis threads to use. Note that in addition to the analysis threads, kanpig keeps one dedicated IO thread
for VCF reading and writing.

# Experimental Parameter Details

These parameters have a varying effect on the results and are not guaranteed to be stable across releases. 

## `--try-exact`
Before performing the path-finding algorithm that applies haplotypes to the variant graph, perform a 1-to-1 comparison
of the haplotypes to each node in the variant graph. If a single node matches above `--sizesim` and `--seqsim`, the
path-finding is skipped and haplotype applied to support the node. 

This parameter will boost the specificity and speed of kanpig at the cost of recall.

## `--prune`
Similar to `--try-exact`, a 1-to-1 comparison is performed before path-finding. If any matches are found, all paths
which do not traverse the matching nodes are pruned from the variant graph. 

This parameter will boost the specificity and speed of kanpig at the cost of recall.

## `--maxhom`

When performing kmer-featurization of sequences (from reads or variants), homopolymer runs above `--maxhom` are trimmed
to `--maxhom`. For example, `--maxhom 5` will only count two four-mers in all homopolymer runs above 5bp.

## `--spanoff`

Don't use this.

kanpig - Kmer ANalysis of PIleups for Genotyping
<img src="https://github.com/ACEnglish/kanpig/raw/develop/imgs/icon.png/" style="width:100px;" align="right" style="vertical-align: middle;"> 
------
A fast tool for genotyping structural variants with long-reads.

# 📥 Install
Binaries are available in [releases](https://github.com/ACEnglish/kanpig/releases).

Alternatively, build from source with:
```
git clone https://github.com/ACEnglish/kanpig
cd kanpig
cargo build --release
# executable in ./target/release/kanpig
```

# 🚀 Quick Start
```
kanpig gt --input variant.vcf.gz --reads alignments.bam --reference ref.fa --out output.vcf
```
See `kanpig -h` for all available parameters, most of which are detailed below.

Kanpig can also create a pileup index from a bam that is smaller and faster for the genotyper to parse. This is useful
for long-term or multiple reanalysis operations like N+1 for a cohort.

```
kanpig plup --bam alignments.bam | bedtools sort | bgzip > alignments.plup.gz
tabix -p bed alignments.plup.gz
```

# ⚠️ Current Limitations
* Kanpig expects sequence resolved SVs. Variants with symbolic alts (e.g. `<DEL>`) and BNDs are not parsed.
* Kanpig only looks at read pileups and does not consider split or soft-clipped alignment information. This means
  variants above ~10kbp should be skipped with the `--sizemax` parameter.

# 🔧 Core Parameter Details

The default parameters are tuned to work generally well for genotyping a single sample's VCF, meaning the variants are
all expected to be present in the sample. For a multi-sample VCF (a.k.a. a project-level VCF), the optimal parameters
are still being determined and will likely be dependent on things such as number of samples in the VCF, merging strategy 
of the variants, and sequencing technology.

### `--bed`
A sorted bed file (`bedtools sort`) that restricts kanpig to only analyzing variants with starts and ends within a single bed entry.

### `--ploidy-bed`
This bed file informs kanpig of special regions within chromosomes that should have non-diploid genotypes. For example, a female
human sample shouldn't have any genotypes on chrY. A male human sample should have hemizygous genotypes on chrY and the
non-pseudoautosomal regions of chrX. The [ploidy_beds/](https://github.com/ACEnglish/kanpig/tree/develop/ploidy_beds) directory 
has example bed files for GRCh38. All regions not within the `--ploidy-bed` (or if no bed is provided) are assumed to be diploid.

### `--chunksize`
Kanpig will build local variant graphs from groups of variants in a 'neighborhood'. These neighborhoods are determined by making the maximum end position
of an upstream neighborhood's variants at least `chunksize` base-pairs away from the next neighborhood's variants' minimum start position.

This chunksize also determines the region over which read pileups are generated. Only reads with at least `mapq` mapping quality, 
passing the `mapflag` filter, and which fully span the neighborhood are considered.

This is an important parameter because too small of a `chunksize` may not recruit distant read pileups which support variants. Similarly, 
too large of a value may create long neighborhoods with many SVs which are also too large for reads to fully-span.

### `--sizemin` and `--sizemax`
Variant sizes are determined by `abs(length(ALT) - length(REF))`. Genotypes of variants not within the size boundaries are set to missing (`./.`).

### `--sizesim` and `--seqsim`
When applying a haplotype to a variant graph, only paths above these two thresholds are allowed. If there are multiple
paths above the threshold, the one with the highest score is kept. Generally, `0.90` is well balanced
whereas lower thresholds will boost recall at the cost of precision and vice versa for higher thresholds. 

### `--gpenalty` and `--fpenalty`
The similarity of a path to the graph is used to to compute a score and the highest score kept. The scoring formula is

```
Score(P) = ((SS + SZ) / 2) − (λg ⋅ ∣L(P)−E∣) - (λf ⋅ N)
``` 

where `SS` and `SZ` are sequence and size similarity,  `L(P)` is the number of nodes in the path, and `E` is the number of 
pileups in the haplotype, and `N` is the number of putative false-negatives in the variant graph. 

The penalty factor `λg` helps reduce paths with split variant representations. The penalty factor `λf` helps penalizes
false-negatives in the variant graph. Details on the scoring penalties are in [the wiki](https://github.com/ACEnglish/kanpig/wiki/Scoring-Function).

### `--maxpaths`
When performing path-finding, this threshold limits the number of paths which are checked. A lower `maxpaths` will
speed up runtime but may come at a cost of recall. A higher `maxpaths` is slower and may come at a cost to
specificity.

### `--hapsim`
After performing kmeans clustering on reads to determine the two haplotypes, if the two haplotypes have a size similarity 
above `hapsim`, they are consolidated into a homozygous allele.

### `--threads`
Number of analysis threads to use. Note that in addition to the analysis threads, kanpig keeps one dedicated IO thread
for VCF reading and writing.

# 📝 Annotations

The `SAMPLE` column fields populated by kanpig are:

| Field   | Description |
|---------|-------------|
| **FT**  | Bit flag for properties of the variant's genotyping. Flags == 0 are considered PASS. |
| **SQ**  | Phred scaled likelihood variant alternate is present in the sample |
| **GQ**  | Phred scale difference between most and second-most likely genotypes |
| **PS**  | Each chunk of variants is assigned a phase set |
| **DP**  | Read coverage over the region |
| **AD**  | Read coverage supporting the reference and alternate alleles. |
| **KS**  | [Kanpig score](https://github.com/ACEnglish/kanpig/wiki/Scoring-Function) |

Details of `FT`
| Flag   | Description |
|--------|-------------|
| 0x1    | The genotype observed from variants matching paths is not equal to the genotype observed from measuring the proportions of reads supporting the two alleles. |
| 0x2    | The genotype quality is less than 5 |
| 0x4    | The depth (DP) is less than 5 |
| 0x8    | The sample quality (SQ) is less than 5 (only present on non-ref variants) |
| 0x16   | The number of reads supporting the alternate allele less than 5 (only present on non-ref variants) |
| 0x32   | The best scoring path through the variant graph only used part of the haplotype. This may be indicative of a false-negative in the variant graph. |

# 🔌 Compute Resources

Kanpig is highly parallelized and will fully utilize all threads it is given. However, hyperthreading doesn't seem to
help and therefore the number of threads should probably be limited to the number of physical processors available. 

For memory, a general rule is kanpig will need about 20x the size of the compressed `.vcf.gz`. The minimum required 
memory is also dependent on the number of threads running as each will need space for its processing. For example, 
a 1.6Gb vcf (~5 million SVs) using 16 cores needs at least 32Gb of RAM. That same vcf with 8 or 4 cores needs at least
 24Gb and 20Gb of RAM, respectively. 

# 🔬 Experimental Parameter Details

These parameters have a varying effect on the results and are not guaranteed to be stable across releases. 

### `--try-exact`
Before performing the path-finding algorithm that applies a haplotype to the variant graph, perform a 1-to-1 comparison
of the haplotype to each node in the variant graph. If a single node matches above `sizesim` and `seqsim`, the
path-finding is skipped and haplotype applied to the node. 

This parameter will boost the specificity and speed of kanpig at the cost of recall.

### `--prune`
Similar to `try-exact`, a 1-to-1 comparison is performed before path-finding. If any matches are found, all paths
which do not traverse the matching nodes are pruned from the variant graph. 

This parameter will boost the specificity and speed of kanpig at the cost of recall.

### `--maxhom`

When performing kmer-featurization of sequences (from reads or variants), homopolymer runs above `maxhom` are trimmed
to `maxhom`. For example, `--maxhom 5` will only count two four-mers in homopolymer runs above 5bp.

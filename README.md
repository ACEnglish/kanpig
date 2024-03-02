kanpig - Kmer ANalysis of PIleups for Genotyping
<img src="https://github.com/ACEnglish/kanpig/raw/develop/imgs/icon.png/" style="width:100px;" align="right" style="vertical-align: middle;"> 
------
A fast tool for genotyping structural variants with long-reads.

Install
-------
```
git clone https://github.com/ACEnglish/kanpig
cd kanpig
cargo build --release
# executable in ./target/release/kanpig
```

Usage
-----
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
      --kmer <KMER>            Kmer size for featurization [default: 4]
      --chunksize <CHUNKSIZE>  Minimum distance between variants to create independent graphs [default: 100]
      --passonly               Only analyze reads with PASS FILTER
      --sizemin <SIZEMIN>      Minimum size of variant to analyze [default: 20]
      --sizemax <SIZEMAX>      Maximum size of variant to analyze [default: 50000]
      --maxpaths <MAXPATHS>    Maximum number of paths in a graph to traverse [default: 1000]
      --seqsim <SEQSIM>        Minimum sequence similarity for paths [default: 0.9]
      --sizesim <SIZESIM>      Minimum size similarity for paths [default: 0.9]
      --minkfreq <MINKFREQ>    Minimum frequency of kmer [default: 1]
  -h, --help                   Print help
  -V, --version                Print version
```

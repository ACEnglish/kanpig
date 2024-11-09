This is a small demo of how to run kanpig.
This genotypes a handful of variants from GIAB v1.1 on grch38 chr5:15903849-16153359 

Files
-----
- run.sh : script to run 
- hg002.test.bam : long-reads
- hg002.test.vcf.gz : example variants
- small.chr5.fa.gz : highly compressed chr5 reference fasta
- summary.sh : optional step to analyze the results
- quick_compare.py : used by summary.sh
- extract.py : script used to make a highly compressed chr5 reference fasta

Setup
-----
The `run.sh` assumes that `../target/release/kanpig` exists. Update the path if the executable is somewhere else.

The `summary.sh` is optional and assumes bcftools, tabix, python3 (with pysam) is available.

Instructions
------------
```bash
bash run.sh
# Optional
bash summary.sh
```

Check `results.vcf` to see the kanpig output. 

The `summary.sh` will merge the original input vcf with the kanpig results and will report how many variants 
were correctly genotyped. Alternatively, open `answer.vcf.gz` to read the genotypes side-by-side.

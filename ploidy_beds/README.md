# ploidy beds

Tab-delimited file with chromosome, start, end and ploidy of regions.
Acceptable ploidys are 0 for regions which shouldn't be genotyped and 1 for hemizygous regions.

The male bed file chrX regions were created using

```bash
bedtools complement -g genome.bed -i par.bed
```

Where `par.bed` is the pseudoautosomal regions as defined 
[here](https://github.com/lh3/dipcall/blob/master/data/hs38.PAR.bed)

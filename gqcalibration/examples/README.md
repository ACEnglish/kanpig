Example GQ Configs
==================
Example gq configs are available and described below. These have not been tested on other sequencing experiments
i.e. using the hifi.38x gq configs on, say, an ONT 20x sequencing experiment. That isn't to say that GQ calibration on
one type of long-read sequencing isn't useful for another type of long-read sequencing, just that these example gq
configs shouldn't be automatic defaults used in your kanpig runs.

GIAB v1.1
---------

Sequencing from <path> was converted to plup via

```
kanpig plup --threads 4 --bam https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam \
    | bedtools sort -header | bgzip > HG002.revio38x.plup.gz
tabix -p bed HG002.revio38x.plup.gz
```

The truth set SV (`stvar`) VCF and high-confidence bed file from GIAB v1.1 was downloaded from
[NIST](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/).

Genotyping and gq calibration performed with

```
vcf=GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz
bed=GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
reads=HG002.revio38x.plup.gz
kanpig gt --input ${vcf} --reads ${reads} \
    --reference GRCh38_1kg_mainchrs.fa --threads 4 \
    | bcftools sort -O z -o kanpig.giab.vcf.gz
tabix kanpig.giab.vcf.gz

bcftools merge -m none --force-samples ${vcf} kanpig.giab.vcf.gz -O z -o giab.merged.vcf.gz
tabix giab.merged.vcf.gz

estimate_params.py giab.merged.vcf.gz giabv1.1.hifi.38x  --bed ${bed} --all --flat-priors --leaveout 0.10
```

HPRC samples
--------------

The GIAB v1.1 HG002 sequencing from above was used as input reads. 

Using the kanpig publication's [HPRC assembly derived SVs](https://zenodo.org/records/14726292), a multi-sample truth set was created for HG002. 
First, all non-HG002 VCFs were consolidated with `bcftools merge -m none`. Second, we removed SVs from the consolidated
VCF with truvari in order to lessen the chances of kanpig applying coverage to a highly similar non-HG002 SV.

```
truvari bench -b hg002.hprc.vcf.gz -c non-hg002.hprc.vcf.gz --pctseq 0.90 --pctsize 0.90 --short --pick multi -o bench/
bcftools merge -m none -0 hg002.hprc.vcf.gz bench/fp.vcf.gz -O u | bcftools view -s HG002 -O z -o hprc.hg002.vcf.gz
tabix hprc.hg002.vcf.gz
```

Kanpig and parameter estimation was run with parameters

```
vcf=hprc.hg002.vcf.gz
bed=GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
reads=HG002.revio38x.plup.gz

kanpig gt --input ${vcf} --reads ${reads} --seqsim 0.85 --maxpaths 1000 --reference GRCh38_1kg_mainchrs.fa --threads 4 \
    | bcftools sort -O z -o kanpig.hprc.vcf.gz
tabix kanpig.hprc.vcf.gz

bcftools merge -m none --force-samples ${vcf} kanpig.hprc.vcf.gz -O z -o hprc.merged.vcf.gz
tabix hprc.merged.vcf.gz

estimate_params.py --bed ${bed} --all --flat-priors --leaveout 0.10 hprc.merged.vcf.gz outputs/hprc.hg002.hifi.38x
```

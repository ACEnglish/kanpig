Example GQ Configs
==================

Example gq configs are available and described below. These have not been tested on other sequencing experiments
i.e. using the hifi.38x gq configs on, say, an ONT 20x sequencing experiment. That isn't to say that GQ calibration on
one type of long-read sequencing isn't useful for another type of long-read sequencing, just that these example gq
configs shouldn't be automatic defaults used in your kanpig runs.

Note that the default genotyper is Beta. In order to calibrate on a different genotyper, an initial kanpig run should
have a default gqconfig provided e.g.

```json
{
  "mode": "Bino",
  "mixture_fractions": [0.33, 0.34, 0.33],
  "means": [0.03, 0.50, 0.97],
  "precisions": [0.0, 0.0, 0.0],
  "calibration_table": [],
  "metadata": {
    "notes": "Manual Default",
  }
}
```

GIAB v1.1
---------

Sequencing from GIAB was converted to plup via

```bash
kanpig plup --threads 4 --bam https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam \
    | bedtools sort -header | bgzip > HG002.revio38x.plup.gz
tabix -p bed HG002.revio38x.plup.gz
```

The truth set SV (`stvar`) VCF and high-confidence bed file from GIAB v1.1 was downloaded from
[NIST](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.019-20241113/).

Kanpig genotyping and parameter estimation was performed via:

```bash
vcf=GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz
bed=GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
reads=HG002.revio38x.plup.gz

kanpig gt --input ${vcf} --reads ${reads} \
    --reference GRCh38_1kg_mainchrs.fa --threads 4 \
    | bcftools sort -O z -o kanpig.giab.vcf.gz
tabix kanpig.giab.vcf.gz

bcftools merge -m none --force-samples ${vcf} kanpig.giab.vcf.gz -O z -o giab.merged.vcf.gz
tabix giab.merged.vcf.gz

estimate_params.py --bed ${bed} --all --flat-priors --leaveout 0.10 giab.merged.vcf.gz giabv1.1.hifi.38x
```

HPRC samples
--------------

The GIAB v1.1 HG002 sequencing from above was used as input reads. 

Using the kanpig publication's [HPRC assembly derived SVs](https://zenodo.org/records/14726292), a multi-sample truth 
set was created for HG002. First, all non-HG002 VCFs were consolidated with `bcftools merge -m none`. Second, we 
removed SVs from the consolidated VCF with truvari in order to lessen the chances of kanpig applying coverage to a 
highly similar non-HG002 SV.

```bash
truvari bench -b hg002.hprc.vcf.gz \
    -c non-hg002.hprc.vcf.gz \
    --pctseq 0.90 \
    --pctsize 0.90 \
    --short \
    --pick multi \
    -o bench/

bcftools merge -m none -0 -O u \
    hg002.hprc.vcf.gz bench/fp.vcf.gz \
    | bcftools view -s HG002 -O z -o hprc.hg002.vcf.gz

tabix hprc.hg002.vcf.gz
```

Kanpig genotyping and parameter estimation was performed via:

```bash
vcf=hprc.hg002.vcf.gz
bed=GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
reads=HG002.revio38x.plup.gz

kanpig gt --input ${vcf} --reads ${reads} \
    --seqsim 0.85 --maxpaths 1000 \
    --reference GRCh38_1kg_mainchrs.fa --threads 4 \
    | bcftools sort -O z -o kanpig.hprc.vcf.gz

tabix kanpig.hprc.vcf.gz

bcftools merge -m none --force-samples -O z -o hprc.merged.vcf.gz ${vcf} kanpig.hprc.vcf.gz 
tabix hprc.merged.vcf.gz

estimate_params.py --bed ${bed} --all --flat-priors --leaveout 0.10 hprc.merged.vcf.gz hprc.hg002.hifi.38x
```

Another HPRC
------------
An example experiment where we use different sequencing experiment and two other HPRC samples is described in
`experiment_hprc.sh`.

Discovery SVs (pending)
----------------------------

* Grab a sniffles discovery VCF from HG002. 
* truvari bench --pick multi, then build a script that will look at the MatchIds to figure out what the true genotype
  should be, place that on the sniffles variant representation. Need to do the pick multi because sniffles may have
  already overmerged non-SV internal differences e.g. GIAB has compound het insertions that only have a SNP difference,
  which sniffles won't (shouldn't) detect and therefore will only have a single variant representation. All FPs are 0/0
* Then do the kanpig/param est.

Platinum Pedigrees (pending)
----------------------------

The truth vcf for NA12878 is good enough. The thing I'm missing here is if/how to get into the trio mode genotyper / GQs
for calibration. At minimum, we can still use germ genotyper without leveraging relatedness to help inform GQs.

Other sequencing wishlist
-------------------------

Different coverages. ONT.


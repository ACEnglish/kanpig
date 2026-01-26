#!/bin/bash

bURL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.020-20250117/

# Include Bed
wget ${bURL}/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed -O include.bed

# Baseline VCF
bcftools view -i "SVLEN >= 50"  -r chr20 -O z -o baseline.vcf.gz ${bURL}/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz
tabix baseline.vcf.gz
# Clean
rm GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz.tbi

# BAM
bam=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam

samtools view -O BAM -o reads.bam ${bam} chr20
samtools index reads.bam
# Clean
rm HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam.bai




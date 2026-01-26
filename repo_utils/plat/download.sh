CHROM=chr20
S3=s3://platinum-pedigree-data

aws s3 --no-sign-request cp s3://platinum-pedigree-data/truthset_v1.2/NA12878_hq_v1.2.svs.bed.gz include.bed.gz
gunzip include.bed.gz
bcftools view \
    -e 'strlen(REF) >= 5 && strlen(ALT) >= 5' \
    -r chr20 \
    s3://platinum-pedigree-data/truthset_v1.2/NA12878_hq_v1.2.svs.vcf.gz \
    -O z -o baseline.vcf.gz
tabix baseline.vcf.gz

#for i in NA12878 NA12891 NA12892
#do
    #bcftools view -r ${CHROM} -O z -o ${i}.vcf.gz ${S3}/variants/assembly-based/dipcall/GRCh38/${i}.dip.vcf.gz
    #tabix ${i}.vcf.gz
    #samtools view -o ${i}.bam -O BAM ${S3}/data/ont/mapped/GRCh38/${i}.minimap2.bam ${CHROM}
    #samtools index ${i}.bam
#done

#bcftools merge -m none -O u NA12878.vcf.gz NA12891.vcf.gz NA12892.vcf.gz \
    #| bcftools norm -N -m-any \
    #| truvari anno svinfo \
    #| bcftools view -i "SVLEN >= 50" -O z -o asm_baseline.vcf.gz
#tabix asm_baseline.vcf.gz

mkdir -p ${OD}/giab_v1.1/
$kanpig \
    gt \
    --sample kanpig \
    --reads ${TESTSRC}/giab_v1.1/reads.bam \
    --input ${TESTSRC}/giab_v1.1/baseline.vcf.gz \
    --threads 4 \
    --reference ${REF} \
    -o ${OD}/giab_v1.1/default_output.vcf

bcftools sort ${OD}/giab_v1.1/default_output.vcf -O z -o ${OD}/giab_v1.1/default_output.vcf.gz
tabix ${OD}/giab_v1.1/default_output.vcf.gz
bcftools merge -m none -O z \
        ${TESTSRC}/giab_v1.1/baseline.vcf.gz \
        ${OD}/giab_v1.1/default_output.vcf.gz \
        -o ${OD}/giab_v1.1/merged.vcf.gz
tabix ${OD}/giab_v1.1/merged.vcf.gz

mkdir -p ${OD}/giab_v1.1/

echo '### Test gt on a BAM ###'
$kanpig \
    gt \
    --sample kanpig \
    --reads ${TESTSRC}/giab_v1.1/reads.bam \
    --input ${TESTSRC}/giab_v1.1/baseline.vcf.gz \
    --threads 4 \
    --reference ${REF} \
    | bcftools sort  -O z -o ${OD}/giab_v1.1/default_output.vcf.gz

tabix ${OD}/giab_v1.1/default_output.vcf.gz
bcftools merge -m none -O z \
        ${TESTSRC}/giab_v1.1/baseline.vcf.gz \
        ${OD}/giab_v1.1/default_output.vcf.gz \
        -o ${OD}/giab_v1.1/default_merged.vcf.gz
tabix ${OD}/giab_v1.1/default_merged.vcf.gz

python ${TESTSRC}/giab_v1.1/calc_perf.py ${OD}/giab_v1.1/default_merged.vcf.gz ${TESTSRC}/giab_v1.1/include.bed

echo '### Test BAM to PLUP ###'
$kanpig \
    plup \
    --bam ${TESTSRC}/giab_v1.1/reads.bam \
    --threads 4 \
    | bedtools sort -header \
    | bgzip > ${OD}/giab_v1.1/HG002.plup.gz
tabix -p bed ${OD}/giab_v1.1/HG002.plup.gz

echo '### Test gt on a PLUP ###'
$kanpig \
    gt \
    --sample kanpig \
    --reads ${OD}/giab_v1.1/HG002.plup.gz \
    --input ${TESTSRC}/giab_v1.1/baseline.vcf.gz \
    --threads 4 \
    --reference ${REF} \
    | bcftools sort  -O z -o ${OD}/giab_v1.1/plup_output.vcf.gz

# I don't like repeating this. Make it a function
tabix ${OD}/giab_v1.1/plup_output.vcf.gz
bcftools merge -m none -O z \
        ${TESTSRC}/giab_v1.1/baseline.vcf.gz \
        ${OD}/giab_v1.1/plup_output.vcf.gz \
        -o ${OD}/giab_v1.1/plup_merged.vcf.gz
tabix ${OD}/giab_v1.1/plup_merged.vcf.gz

python ${TESTSRC}/giab_v1.1/calc_perf.py ${OD}/giab_v1.1/plup_merged.vcf.gz ${TESTSRC}/giab_v1.1/include.bed

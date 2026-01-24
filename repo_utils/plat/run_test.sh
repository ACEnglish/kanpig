mkdir -p ${OD}/plat/

echo '### Test trios ###'

$kanpig \
    trio \
    --bed ${TESTSRC}/plat/include.bed \
    --proband ${TESTSRC}/plat/NA12878.bam \
    --mother ${TESTSRC}/plat/NA12892.bam \
    --father ${TESTSRC}/plat/NA12891.bam \
    --input ${TESTSRC}/plat/baseline.vcf.gz \
    --threads 4 \
    --reference ${REF} \
    | bcftools sort  -O z -o ${OD}/plat/default_output.vcf.gz

tabix ${OD}/plat/default_output.vcf.gz

bcftools merge -m none -O z \
        ${TESTSRC}/plat/baseline.vcf.gz \
        ${OD}/plat/default_output.vcf.gz \
        -o ${OD}/plat/default_merged.vcf.gz
tabix ${OD}/plat/default_merged.vcf.gz

python ${TESTSRC}/calc_trio_perf.py ${OD}/plat/default_merged.vcf.gz ${TESTSRC}/giab_v1.1/include.bed

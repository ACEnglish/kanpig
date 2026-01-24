mkdir -p ${OD}/mims/

echo '### Test trios ###'

$kanpig \
    mosaic \
    --bed ${TESTSRC}/mims/include.bed \
    --reads ${TESTSRC}/mims/HapMap.cram \
    --input ${TESTSRC}/mims/baseline.vcf.gz \
    --threads 4 \
    --reference ${REF} \
    | bcftools sort  -O z -o ${OD}/mims/default_output.vcf.gz

tabix ${OD}/mims/default_output.vcf.gz

bcftools merge -m none -O z \
        ${TESTSRC}/mims/baseline.vcf.gz \
        ${OD}/mims/default_output.vcf.gz \
        -o ${OD}/mims/default_merged.vcf.gz
tabix ${OD}/mims/default_merged.vcf.gz

python ${TESTSRC}/mims/calc_perf.py ${OD}/mims/default_merged.vcf.gz ${TESTSRC}/mims/include.bed

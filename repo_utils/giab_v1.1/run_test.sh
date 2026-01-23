# Remember that these tests are run from the repo root directory
mkdir -p ${OD}/giab_v1.1/
run giab_default cargo run -- \
    gt \
    --sample kanpig \
    --reads ${TESTSRC}/giab_v1.1/reads.bam \
    --input ${TESTSRC}/giab_v1.1/baseline.vcf.gz \
    --threads 4 \
    --reference ${REF} \
    -o ${OD}/giab_v1.1/default_output.vcf


if [ $giab_default ]; then
    assert_exit_code 0

    bcftools sort ${OD}/giab_v1.1/default_output.vcf -O z -o ${OD}/giab_v1.1/default_output.vcf.gz
    tabix ${OD}/giab_v1.1/default_output.vcf.gz
    bcftools merge -m none -O z \
            ${TESTSRC}/giab_v1.1/baseline.vcf.gz \
            ${OD}/giab_v1.1/default_output.vcf.gz \
            -o ${OD}/giab_v1.1/merged.vcf.gz
    tabix ${OD}/giab_v1.1/merged.vcf.gz

    # This needs some asserts
    # I wouldn't mind keeping commit to perf records, either
    # Wrap this in an actual test
    python ${TESTSRC}/giab_v1.1/calc_perf.py ${OD}/giab_v1.1/merged.vcf.gz ${TESTSRC}/giab_v1.1/include.bed
    assert_exit_code 0
fi

# Be sure to exit back out
cd - 

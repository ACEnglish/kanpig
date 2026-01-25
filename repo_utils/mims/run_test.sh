OD=$1
SRC=$2

echo '### Test mosaic ###'

$kanpig \
    mosaic \
    --bed ${SRC}/include.bed \
    --reads ${SRC}/SMAFIKCF6M5Z.cram \
    --sample samp0 \
    --reads ${SRC}/SMAFIOYCBUCR.cram \
    --sample samp1 \
    --reads ${SRC}/SMAFIR3J8UVK.cram \
    --sample samp2 \
    --input ${SRC}/baseline.vcf.gz \
    --threads 4 \
    --reference ${REF} \
    | bcftools sort  -O z -o ${OD}/default_output.vcf.gz

tabix ${OD}/default_output.vcf.gz

bcftools merge -m none -O z \
        ${SRC}/baseline.vcf.gz \
        ${OD}/default_output.vcf.gz \
        -o ${OD}/default_merged.vcf.gz
tabix ${OD}/default_merged.vcf.gz

python ${SRC}/calc_perf.py ${OD}/default_merged.vcf.gz ${SRC}/include.bed

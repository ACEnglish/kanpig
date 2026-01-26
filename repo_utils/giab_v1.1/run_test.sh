OD=$1
SRC=$2

echo '### Test gt on a BAM ###'
$kanpig \
    gt \
    --sample kanpig \
    --reads ${SRC}/reads.bam \
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

echo '### Test BAM to PLUP ###'
$kanpig \
    plup \
    --bam ${SRC}/reads.bam \
    --threads 4 \
    | bedtools sort -header \
    | bgzip > ${OD}/HG002.plup.gz
tabix -p bed ${OD}/HG002.plup.gz

echo '### Test gt on a PLUP ###'
$kanpig \
    gt \
    --sample kanpig \
    --reads ${OD}/HG002.plup.gz \
    --input ${SRC}/baseline.vcf.gz \
    --threads 4 \
    --reference ${REF} \
    | bcftools sort  -O z -o ${OD}/plup_output.vcf.gz

# I don't like repeating this. Make it a function
tabix ${OD}/plup_output.vcf.gz
bcftools merge -m none -O z \
        ${SRC}/baseline.vcf.gz \
        ${OD}/plup_output.vcf.gz \
        -o ${OD}/plup_merged.vcf.gz
tabix ${OD}/plup_merged.vcf.gz

python ${SRC}/calc_perf.py ${OD}/plup_merged.vcf.gz ${SRC}/include.bed

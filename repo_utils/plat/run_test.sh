OD=$1
SRC=$2

echo '### Test trio ###'

$kanpig \
    trio \
    --bed ${SRC}/include.bed \
    --proband ${SRC}/NA12878.bam \
    --mother ${SRC}/NA12892.bam \
    --father ${SRC}/NA12891.bam \
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

bcftools +mendelian2 ${OD}/default_output.vcf.gz -p 1X:PRO,PAT,MAT -m c

OD=$1
SRC=$2

$kanpig gt \
    --reads ${SRC}/HG002_withHP.plup.gz \
    --input ${SRC}/baseline.vcf.gz \
    --reference ${REF} \
    -o ${OD}/with.vcf

result=$(bcftools query -f "[%PS]\n" ${OD}/with.vcf)
test "$result" = "8903002" || exit 1


$kanpig gt \
    --reads ${SRC}/HG002_woHP.plup.gz \
    --input ${SRC}/baseline.vcf.gz \
    --reference ${REF} \
    -o ${OD}/without.vcf

result=$(bcftools query -f "[%PS]\n" ${OD}/without.vcf)
test "$result" = "9015824" || exit 1

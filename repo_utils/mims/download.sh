# I need to get the full coverage
# Grabe a few other BCM samples so that you can also test the multi-sample parsing
for fn in f346499b-3f45-4f3b-a2ac-7d15351fe6e2/@@download/SMAFIOYCBUCR.cram \
          5b8b7376-00b9-4621-bcd8-f36ebdca4cfc/@@download/SMAFIKCF6M5Z.cram \
          3ed3f0a7-0a98-4f38-bd5e-219103818416/@@download/SMAFIR3J8UVK.cram
do
    name=$(basename ${fn})
    samtools view -O CRAM -o ${name} https://${ACCESS_KEY_ID}:${SECRET_ACCESS_KEY}@data.smaht.org/output-files/${fn} chr20
    samtools index ${name}
done

bcftools view -r chr20 https://github.com/BCM-HGSC/SMaHT_MIMS/raw/refs/heads/main/benchmark_v2/smaht_mims_sv_v2_easy.vcf.gz -O z -o baseline.vcf.gz
tabix baseline.vcf.gz
rm smaht_mims_sv_v2_easy.vcf.gz.tbi

wget -O include.bed https://github.com/BCM-HGSC/SMaHT_MIMS/raw/refs/heads/main/benchmark_v2/smaht_mims_sv_v2_easy_regions.bed

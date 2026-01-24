samtools view -O CRAM -o HapMap.cram https://${ACCESS_KEY_ID}:${SECRET_ACCESS_KEY}@data.smaht.org/output-files/5b8b7376-00b9-4621-bcd8-f36ebdca4cfc/@@download/SMAFIKCF6M5Z.cram	chr20
samtools index HapMap.cram
rm SMAFIKCF6M5Z.cram.crai

bcftools view -r chr20 https://github.com/BCM-HGSC/SMaHT_MIMS/raw/refs/heads/main/benchmark_v2/smaht_mims_sv_v2_easy.vcf.gz -O z -o baseline.vcf.gz
tabix baseline.vcf.gz
rm smaht_mims_sv_v2_easy.vcf.gz.tbi

wget -O include.bed https://github.com/BCM-HGSC/SMaHT_MIMS/raw/refs/heads/main/benchmark_v2/smaht_mims_sv_v2_easy_regions.bed

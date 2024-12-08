bcftools sort result.vcf -O z -o result.vcf.gz
tabix result.vcf.gz
bcftools merge hg002.test.vcf.gz result.vcf.gz -O z -o answer.vcf.gz
tabix answer.vcf.gz
python quick_compare.py

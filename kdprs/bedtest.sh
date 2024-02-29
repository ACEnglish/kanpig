#./target/release/kdprs \
cargo run --release -- \
    --input test/test2.vcf.gz \
    --bam /Users/english/code/kfdphase/test/NA24385.chr20.bam \
    --reference /Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa \
    --bed test/test.chr20.bed \
    -o test/nope.vcf

cut -f1-10 test/nope.vcf | bcftools sort -O z -o  test/nope.vcf.gz
tabix test/nope.vcf.gz
rm -rf test/new
truvari bench --includebed test/test.chr20.bed -b ../test/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz -c test/nope.vcf.gz -o test/new

rm -rf test/new_noref
truvari bench --includebed test/test.chr20.bed -b ../test/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz -c test/nope.vcf.gz --no-ref a -o test/new_noref

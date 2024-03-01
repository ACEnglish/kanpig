set -e
#./target/release/kdprs \
cargo run --release -- \
    --input test/test2.vcf.gz \
    --bam /Users/english/code/kfdphase/test/NA24385.chr20.bam \
    --reference /Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa \
    --bed /Users/english/code/kfdphase/test/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed \
    --sizesim 0.95 --seqsim 0.95 --threads 4 \
    -o test/fast.vcf
# --bam test/GIABHG002.bam \

#cut -f1-10 test/hc.vcf | bcftools sort -O z -o  test/hc.vcf.gz
#tabix test/hc.vcf.gz
#truvari bench --includebed /Users/english/code/kfdphase/test/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed \
        #-b ../test/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
        #-c test/hc.vcf.gz -o test/hcbench_all/

#truvari bench --includebed /Users/english/code/kfdphase/test/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed \
        #-b ../test/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
        #-c test/hc.vcf.gz --no-ref a -o test/hcbench_all/

# truvari refine -f ~/code/references/grch38/GRCh38_1kg_mainchrs.fa -U -u -R --regions new_noref/candidate.refine.bed new_noref

set -e
bed=test_rs/test.chr20.bed

create() {
    #./target/release/kdprs \
    cargo run --release -- \
        --input test_rs/test2.vcf.gz \
        --bam /Users/english/code/kanplug/experiments/test_rs/NA24385.chr20.bam \
        --reference /Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa \
        --bed $bed \
        --sizesim 0.90 --seqsim 0.90 --threads 4 \
        --maxpaths 5000 \
        -o test_rs/hc.vcf
    #--bed /Users/english/code/kfdphase/test/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed \
    # --bam /Users/english/code/kanplug/experiments/test_rs/GIABHG002.bam \
}

bench() {
    cut -f1-10 test_rs/hc.vcf | bcftools sort -O z -o  test_rs/hc.vcf.gz
    tabix test_rs/hc.vcf.gz
    rm -rf test_rs/hcbench_all
    truvari bench --includebed $bed \
        -b test_rs/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
        -c test_rs/hc.vcf.gz -o test_rs/hcbench_all/

    rm -rf test_rs/hcbench_noref
    truvari bench --includebed $bed \
        -b test_rs/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
        -c test_rs/hc.vcf.gz --no-ref a -o test_rs/hcbench_noref/

    truvari refine -f ~/code/references/grch38/GRCh38_1kg_mainchrs.fa -U -u -R \
            --regions test_rs/hcbench_noref/candidate.refine.bed test_rs/hcbench_noref
}

create

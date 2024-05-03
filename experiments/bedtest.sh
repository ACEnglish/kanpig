set -e
bed=test_rs/test.chr20.bed

create() {
    #../target/release/kanpig \
    time cargo run --release -- \
        --input test_rs/test2.vcf.gz \
        --bam /Users/english/code/kanpig/experiments/test_rs/NA24385.chr20.bam \
        --reference /Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa \
        --sizemin 20 \
        --sizesim 0.95 --seqsim 0.90 --threads 4 \
        --maxpaths 10000 --mapq 20 --hapsim 0.98 \
        --chunksize 1000 --prune --try-exact \
        -o test_rs/hc.vcf --bed $bed 
    # --bed /Users/english/code/kanpig/test/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed \
    # --bam /Users/english/code/kanpig/experiments/test_rs/GIABHG002.bam \
}

bench_lite() {
    rm -rf test_rs/hcbench_all
    truvari bench --includebed $bed \
        -b test_rs/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
        -c test_rs/hc.vcf.gz -o test_rs/hcbench_all/ \
        --pctsize 0.90 --pctseq 0.90
}

bench_medium() {
    rm -rf test_rs/hcbench_noref
    truvari bench --includebed $bed \
        -b test_rs/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
        -c test_rs/hc.vcf.gz --no-ref a -o test_rs/hcbench_noref/ \
        --pctsize 0.90 --pctseq 0.90
}

bench_full() {
    truvari refine -f ~/code/references/grch38/GRCh38_1kg_mainchrs.fa -U -u -R \
            --regions test_rs/hcbench_noref/candidate.refine.bed test_rs/hcbench_noref
}

create
bcftools sort -O z -o test_rs/hc.vcf.gz test_rs/hc.vcf
tabix test_rs/hc.vcf.gz
#bench_lite
bench_medium
#bench_full

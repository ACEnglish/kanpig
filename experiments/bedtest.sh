set -e
bed=test_rs/test.chr20.bed

create() {
    #time ../target/release/kanpig \
    #time kanpig-v0.2.0-x86_64-apple-darwin/kanpig \
    time cargo run --release -- \
        --input test_rs/chr20.vcf.gz \
        --bam /Users/english/code/kanpig/experiments/test_rs/NA24385.chr20.bam \
        --reference /Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa \
        --sizemin 50 \
        --sizesim 0.95 --seqsim 0.90 --threads 5 \
        --maxpaths 1000 --mapq 20 --hapsim 0.98 \
        --chunksize 100 --maxhom 0 \
        --sample doesthiswork \
        --bed $bed -o test_rs/hc.vcf
            #| bcftools sort -O z -o test_rs/hc.vcf.gz
    # --bed /Users/english/code/kanpig/test/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed \
    # --bam /Users/english/code/kanpig/experiments/test_rs/GIABHG002.bam \
}

create_zero() {
    #time ../target/release/kanpig \
    #time kanpig-v0.2.0-x86_64-apple-darwin/kanpig \
    time cargo run --release -- \
        --input test_rs/chr20.vcf.gz \
        --bam /Users/english/code/kanpig/experiments/test_rs/NA24385.chr20.bam \
        --reference /Users/english/code/references/grch38/GRCh38_1kg_mainchrs.fa \
        --sizemin 50 \
        --sizesim 0.95 --seqsim 0.90 --threads 5 \
        --maxpaths 1000 --mapq 20 --hapsim 0.98 \
        --chunksize 100 --maxhom 0 \
        --sample doesthiswork \
        --factor 0.03 \
        --bed $bed -o test_rs/hc_zero.vcf
            #| bcftools sort -O z -o test_rs/hc.vcf.gz
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
    oname=${1}_bench/
    rm -rf $oname
    truvari bench --includebed $bed \
        -b test_rs/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
        -c $1 --no-ref a -o $oname \
        --pctsize 0.90 --pctseq 0.90 --pick ac
}

bench_full() {
    oname=${1}_bench/
    truvari refine -f ~/code/references/grch38/GRCh38_1kg_mainchrs.fa -U -u -R \
            --regions $oname/candidate.refine.bed $oname
}

create
bcftools sort -O z -o test_rs/hc.vcf.gz test_rs/hc.vcf
tabix test_rs/hc.vcf.gz
bench_medium test_rs/hc.vcf.gz
#bench_full test_rs/hc.vcf.gz

#create_zero
#bcftools sort -O z -o test_rs/hc_zero.vcf.gz test_rs/hc_zero.vcf
#tabix test_rs/hc_zero.vcf.gz
#bench_medium  test_rs/hc_zero.vcf.gz
#bench_full test_rs/hc_zero.vcf.gz

pvcf=non-hg002.hprc.vcf.gz
reference=GRCh38_1kg_mainchrs.fa
bed_file=GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed

set -e

extract_collapse() {
    sample=$1
    bcftools view -s ${1} ${pvcf} -O z -o ${sample}.vcf.gz
    tabix ${sample}.vcf.gz

    truvari collapse --dynthresh 5,50,50,1500 --keep common -i ${sample}.vcf.gz -c removed.vcf.gz --gt all \
        | bcftools +setGT - -- -t . -n 0p \
        | bcftools sort -O z -o ${sample}.collapsed.vcf.gz
    tabix ${sample}.collapsed.vcf.gz

    
}

train_sample=HG01358
test_sample=HG01978
extract_collapse ${train_sample}
extract_collapse ${test_sample}

# Run once on default with the first sample
kanpig gt --input ${train_sample}.collapsed.vcf.gz --reads ${train_sample}.plup.gz --seqsim 0.85 --threads 4  \
    --reference ${reference} \
    | bcftools sort -O z -o kanpig.vcf.gz
tabix kanpig.vcf.gz
bcftools merge -m none --force-samples -O z -o ${train_sample}.merged.vcf.gz ${train_sample}.collapsed.vcf.gz kanpig.vcf.gz 
tabix ${train_sample}.merged.vcf.gz

# Calibrate
python estimate_params.py ${train_sample}.merged.vcf.gz hprcTrain \
    --all --leaveout 0.10 --bed ${bed_file}

# Run the second sample once with default gqconfig
kanpig gt --input ${test_sample}.collapsed.vcf.gz --reads ${test_sample}.plup.gz --seqsim 0.85 --threads 4  \
    --reference ${reference} \
    | bcftools sort -O z -o kanpig.vcf.gz
tabix kanpig.vcf.gz
bcftools merge -m none --force-samples -O z -o ${test_sample}.merged.vcf.gz ${test_sample}.collapsed.vcf.gz kanpig.vcf.gz 
tabix ${test_sample}.merged.vcf.gz

# Grab its stats
python estimate_params.py ${train_sample}.merged.vcf.gz hprcTest_default \
    --all --no-calibrate --bed ${bed_file}

# Run the second sample again with test_sample gqconfig
kanpig gt --input ${test_sample}.collapsed.vcf.gz --reads ${test_sample}.plup.gz --seqsim 0.85 --threads 4  \
    --reference ${reference} --gqconfig hprcTrain_gqconfig.json \
    | bcftools sort -O z -o kanpig.vcf.gz
tabix kanpig.vcf.gz
bcftools merge -m none --force-samples -O z -o ${test_sample}.gq.merged.vcf.gz ${test_sample}.collapsed.vcf.gz kanpig.vcf.gz 
tabix ${test_sample}.gq.merged.vcf.gz

python estimate_params.py ${test_sample}.gq.merged.vcf.gz hprcTest_withgq \
    --all --no-calibrate --bed ${bed_file}

# Run the second sample again with the correct calibration
# The idea being that this should have pretty bad GQs because the calibration isn't for this depth
kanpig gt --input ${test_sample}.collapsed.vcf.gz --reads ${test_sample}.plup.gz --seqsim 0.85 --threads 4  \
    --reference ${reference} --gqconfig hprcTrain_gqconfig.json \
    | bcftools sort -O z -o kanpig.vcf.gz
tabix kanpig.vcf.gz
bcftools merge -m none --force-samples -O z -o ${test_sample}.gq.merged.vcf.gz ${test_sample}.collapsed.vcf.gz kanpig.vcf.gz 
tabix ${test_sample}.gq.merged.vcf.gz

python estimate_params.py ${test_sample}.gq.merged.vcf.gz hprcTest_withgq \
    --summary --bed ${bed_file} 

# Let's see what happens when we use the same hprc HG002 model which used higher coverage data
kanpig gt --input ${test_sample}.collapsed.vcf.gz --reads ${test_sample}.plup.gz --seqsim 0.85 --threads 4  \
    --reference ${reference} --gqconfig hprc.hg002.hifi.38x_gqconfig.json \
    | bcftools sort -O z -o kanpig.vcf.gz
tabix kanpig.vcf.gz
bcftools merge -m none --force-samples -O z -o ${test_sample}.badgq.merged.vcf.gz ${test_sample}.collapsed.vcf.gz kanpig.vcf.gz 
tabix ${test_sample}.badgq.merged.vcf.gz

python estimate_params.py ${test_sample}.badgq.merged.vcf.gz hprcTest_withbadgq \
    --summary --bed GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed 

# And, lets see what happens with a very different, single-sample model
kanpig gt --input ${test_sample}.collapsed.vcf.gz --reads ${test_sample}.plup.gz --seqsim 0.85 --threads 4  \
    --reference ${reference} \
    --gqconfig giabv1.1.hifi.38x_gqconfig.json \
    | bcftools sort -O z -o kanpig.vcf.gz
tabix kanpig.vcf.gz
bcftools merge -m none --force-samples -O z -o ${test_sample}.realbadgq.merged.vcf.gz ${test_sample}.collapsed.vcf.gz kanpig.vcf.gz 
tabix ${test_sample}.realbadgq.merged.vcf.gz

python estimate_params.py ${test_sample}.realbadgq.merged.vcf.gz hprcTest_withrealbadgq \
    --summary --bed ${bed_file} 



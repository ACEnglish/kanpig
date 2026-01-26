

# GRCh38 chr20 reference
samtools faidx https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz \
    chr20 > ref.fa
samtools faidx ref.fa

# GIAB v1.1

bash giab_v1.1/download.sh

# Platinum Pedigres

bash plat/download.sh

# MIMS

bash mims/download.sh



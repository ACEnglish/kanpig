#def extract_large_insertions(reads):
chrom, start, end = "chr20", 20827970, 20827980 # insertion
all_cov = 0
for pileup_column in bam.pileup(chrom, start, end, truncate=True):
    # Check for deletions
    all_cov += pileup_column.n
    for pileup_read in pileup_column.pileups:
        if pileup_read.indel > 20:  # Insertion greater than 20 bp
            print(pileup_column.reference_pos, pileup_read.indel)
            print(pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+pileup_read.indel])
print(all_cov / (end - start))

#def extract_large_deletions(reads):
chrom, start, end = "chr20", 28797607-50, 28797607+120+50 # deletion
all_cov = 0
for pileup_column in bam.pileup(chrom, start, end, truncate=True):
    # Check for deletions
    all_cov += pileup_column.n
    for pileup_read in pileup_column.pileups:
        if pileup_read.indel < -20:  # Deletion greater than 20 bp
            print(pileup_column.reference_pos, pileup_read.indel)
print(all_cov / (end - start))

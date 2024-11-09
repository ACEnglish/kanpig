import pysam

fasta = pysam.FastaFile("chr5.fa.gz")
start = 13903849 - 100
end = 14153359 + 100
seq = fasta.fetch('chr5', start, end)
print(">chr5")
print("N" * start, seq, "N" * (181538259 - end))

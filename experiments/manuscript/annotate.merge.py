"""
Given a genotyped VCF that has been bcftools merged to its input high-confidence VCF,
annotate each VCF entry with if the genotype is correct
"""
import sys
import pysam
import truvari

# I should be subsetting to the confident bed.

v = pysam.VariantFile(sys.argv[1])
header = v.header.copy()
header.add_line(('##INFO=<ID=GTState,Number=1,Type=String,'
                 'Description="Full genotype comparison state">'))
header.add_line(('##INFO=<ID=ALState,Number=.,Type=String,'
                 'Description="Each Allele State">'))
o = pysam.VariantFile("/dev/stdout", 'w', header=header)
for entry in v:
    if truvari.entry_size(entry) > 10000:
        continue
    entry.translate(header)
    if None in entry.samples[0]["GT"]:
        entry.info["GTState"] = "Missing"
        o.write(entry)
        continue
    if None in entry.samples[1]["GT"]:
        entry.info["GTState"] = "Filtered"
        o.write(entry)
        continue
    c = sorted(list(entry.samples[0]['GT']))
    b = sorted(list(entry.samples[1]['GT']))
    entry.info["GTState"] = "Concordant" if b == c else "Discordant"
    astate = []
    for i,j in zip(b, c):
        state = "T" if i == j else "F"
        condition = "P" if j == 1 else "N"
        astate.append(state + condition)
    entry.info["ALState"] = astate
    o.write(entry)
    


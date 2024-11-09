import pysam

vcf = pysam.VariantFile("answer.vcf.gz")

correct = 0
count = 0
for entry in vcf:
    if None in entry.samples[1]['GT'] or None in entry.samples[0]['GT']:
        continue
    bcnt = sum([_ if _ else 0 for _ in entry.samples[0]['GT']])
    ccnt = sum([_ if _ else 0 for _ in entry.samples[1]['GT']])
    count += 1
    if bcnt == ccnt:
        correct += 1
print('result:', correct, '/', count, '=', correct/count*100, '%')

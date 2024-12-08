import pysam

vcf = pysam.VariantFile("answer.vcf.gz")

correct = 0
count = 0
for entry in vcf:
    bcnt = sum(entry.samples[0]['GT'])
    ccnt = sum(entry.samples[1]['GT'])
    count += 1
    if bcnt == ccnt:
        correct += 1

print('result:', correct, '/', count, '=', correct/count*100, '%')

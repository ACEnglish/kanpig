OD=$1
SRC=$2

echo '### Test gqcalib/bindings'

maturin develop --release --features python

python gqcalibration/estimate_params.py ${SRC}/variants.vcf.gz ${OD}/calib

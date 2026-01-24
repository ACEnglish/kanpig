
mkdir -p ${OD}/pybind
echo '### Test gqcalib/bindings'

maturin develop --release --features python

python gqcalibration/estimate_params.py ${TESTSRC}/pybind/variants.vcf.gz ${OD}/calib

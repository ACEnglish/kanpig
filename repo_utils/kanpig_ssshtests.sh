# Functional tests

# Run from truvari's base directory
cd "$( dirname "${BASH_SOURCE[0]}" )"/../

GITHASH=$(git rev-parse --short HEAD)
DATE=$(date +%Y%m%d_%H%M%S)
echo "### Testing kanpig commit ${GITHASH} on ${DATE}"

OD=test_results
TESTSRC=repo_utils/
REF=$TESTSRC/GRCh38_chr20.fa
kanpig="cargo run --release -- "

LOGDIR=${TESTSRC}/history/${DATE}_${GITHASH}
mkdir -p ${LOGDIR}
exec > >(tee -a ${LOGDIR}/main.out)
exec 2> >(tee -a ${LOGDIR}/main.err >&2)

# Reset test results
rm -rf $OD
mkdir -p $OD

source $TESTSRC/giab_v1.1/run_test.sh
source $TESTSRC/plat/run_test.sh
source $TESTSRC/mims/run_test.sh
source $TESTSRC/pybind/run_test.sh

# Functional tests
test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
source ssshtest

# Run from truvari's base directory
cd "$( dirname "${BASH_SOURCE[0]}" )"/../

OD=test_results
TESTSRC=repo_utils/
REF=$TESTSRC/GRCh38_chr20.fa

kanpig="cargo run -- "

# Reset test results
rm -rf $OD
mkdir -p $OD

source $TESTSRC/giab_v1.1/run_test.sh
#source $TESTSRC/plat/run_test.sh
#source $TESTSRC/mims/run_test.sh

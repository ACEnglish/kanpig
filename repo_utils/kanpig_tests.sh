# Functional tests
set -e

TESTS=(
    "giab_v1.1"
    "plat"
    "mims"
    "pybind"
)

# Run from truvari's base directory
cd "$( dirname "${BASH_SOURCE[0]}" )"/../

GITHASH=$(git rev-parse --short HEAD)
DATE=$(date +%Y%m%d_%H%M%S)
echo "### Testing kanpig commit ${GITHASH} on ${DATE}"

# Variables needed by sub-tests
export TESTSRC=repo_utils/
export REF=$TESTSRC/GRCh38_chr20.fa
export kanpig="cargo run --release -- "

LOGDIR=${TESTSRC}/history/${DATE}_${GITHASH}
mkdir -p ${LOGDIR}
exec > >(tee -a ${LOGDIR}/main.out)
exec 2> >(tee -a ${LOGDIR}/main.err >&2)

# Function to run a single test
run_test() {
    local test_name=$1
    local test_script="$TESTSRC/$test_name/run_test.sh"
    local out_dir=$LOGDIR/$test_name

    # Reset test results for re-runs
    rm -rf $out_dir
    mkdir -p $out_dir

    if [[ -f "$test_script" ]]; then
        echo "### Running test: $test_name"
        bash $test_script $out_dir $TESTSRC/$test_name || {
            echo "!!! Error: Test '$test_name' failed with exit code $?"
            exit 1
        }
    else
        echo "!!! Error: Test script not found: $test_script"
        return 1
    fi
}

# Main logic
if [[ $# -eq 0 ]]; then
    # No parameters - run all tests
    for test in "${TESTS[@]}"; do
        run_test "$test"
    done
else
    # Parameter given - run specific test
    test_name=$1

    # Validate test name
    if [[ ${TESTS[@]} =~ ${test_name} ]]; then
        run_test "$test_name"
    else
        echo "!!! Error: Unknown test '$test_name'"
        echo "!!! Available tests: ${TESTS[*]}"
        exit 1
        break
    fi
fi

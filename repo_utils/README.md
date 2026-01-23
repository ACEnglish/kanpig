Functional Testing Infrastructure


Setup
=====

Download the baseline variants with `bash download.sh`.
This will take single chromosomes from different benchmarks.

Running
=======

Use `bash kanpig_tests.sh`. 

Details
=======

Each sub-directory focuses on a particular test. They'll have a `test.sh` file that can be run.
Some tests will use remote BAMs/CRAMs but can be setup...
Time can be saved by pre-downloading the relevant test files. This can be achieved by running `bash download.sh`.


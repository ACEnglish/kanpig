v0.1.2
*May 5, 2024*

* New optional hompolymer filter doesn't kmerize long homopolymers
* Improved logging info
* Correcting GQ field
* Correcting kmer counting
* Small speed/memory/io improvements
  * Off-loaded annotation work from the single writer thread to the worker threads and using a large
  multiple of page size for the BufWriter capacity
  * Fewer bam file opens

v0.1.1
*Apr 11, 2024*

* The `--no-prune` flag has been changed to `--prune` since not pruning is a better default.
* Partial haplotypes now only allow up to 3 false negatives for regions with fewer than 500 pileups. More than 500 do
  not attempt partials.
* Partial haplotypes now respect the `--kmer` option.

v0.1.0
*Apr 9, 2024*

Initial version. Works well enough to freeze.

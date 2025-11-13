GQ Calibration
==============


Kanpig uses a beta binomial distribution to build genotype quality scores. The defaults will work well for a
set of SVs derived from a single sample's assembly with higher sequencing coverage.

However, different experiments (e.g. single sample discovery, multi-sample merged) may need their GQs calibrated.

To build a custom calibration, first a truth-set VCF must be acquired and run through kanpig without a `--gqconfig`.
Next, the truth-set and kanpig output VCFs should be merged with 

```
bcftools merge -m none truth.vcf.gz kanpig.vcf.gz -O z -o merged.vcf.gz
```
From this VCF, we create a csv of information needed for the calibration via:

```
python make_param_df.py merged.vcf.gz data.csv
```

Next, the kanpig genotyper must be available to your python environment by building the python bindings
```
maturin develop --release --features python
```

Finally, the gqconfig can be created via
```
python estimate_params.py data.csv config.json
```

This config can now be passed into `kanpig gt --gqconfig config.json`

Example Configs
===============

Two example configs are available. 

The first, `gq.giab.json`, was built off of GIAB v1.1 

The second, `gq.hprc.json`, was built off of just HG002 and zenodo kanipig paper assemblies. We took just HG002
genotypes, intersected those with the merge of all the rest of the samples, and then consolidated the HG002 with the
fn.vcf.gz. This allowed us to preserve the genotype quality scores but spike-in reference homozygous variants.
Collapsing the variants would have been an option, but only if we could have done extra work to keep the baseline GT in
there.

TODOs
=====
* Paths to GIAB and kanpig example data
* Clean up all this code and documentation
  * Hard paths
* Pull in the plots from the notebook to automatically generate those
* Generate GT accuracy reports, too. There's patterns to what should be expected from the raw GTs' accuracy that could
  be useful to generate warnings e.g. <75% genotyping accuracy is a problem.
* Expand on what it means to make a truth-set (have to think about merging)

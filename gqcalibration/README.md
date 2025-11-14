GQ Calibration
==============

Kanpig uses a beta binomial distribution to build genotype quality scores. The defaults will work well for a
set of SVs derived from a single sample's assembly with higher sequencing coverage.

However, different experiments (e.g. single sample discovery, multi-sample merged) may need their GQs calibrated.

The script `estimate_params.py` uses a maximum likelihood estimation on a set of genotypes to set the model's parameters
as well as calibrate genotype quality scores and generate informative plots.

To start, the kanpig python-bindings must be installed by running from this repository's root directory the command

```
maturin develop --release --features python
```

Next, a truth-set VCF must be acquired and run through kanpig on a sequencing experiment (withou a `--gqconfig`).
This is the hardest part of GQ calibration. Example VCFs and how they were derived are in the `examples/` directory.

Once the VCF has been genotyped, we reunite the kanpig results with the truth-set genotypes via:

```
bcftools merge -m none truth.vcf.gz kanpig.vcf.gz -O z -o merged.vcf.gz
```
Finally, we generate the gqconfig with the command:

```
python estimate_params.py merged.vcf.gz my_config
```

This config can now be passed into `kanpig gt --gqconfig my_config.json`

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
* Paths to example data

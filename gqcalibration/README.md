GQ Calibration
==============

Kanpig uses a beta binomial distribution to build genotype quality (GQ) scores. The default parameters for the model 
will produce good genotypes, but may not result in informative GQs. Therefore, different experiments (e.g. single 
sample discovery, multi-sample merged, varying sequencing coverage) may need their GQs calibrated.

The script `estimate_params.py` uses a maximum likelihood estimation on a set of genotypes to set the model's parameters
as well as calibrate genotype quality scores and generate informative plots.

To start, kanpig's python bindings must be installed by using maturin from this repository's root directory

```bash
maturin develop --release --features python
```

Other dependencies for the python script are in `requirements.txt` and can be installed with `pip`.

A truth-set VCF must be acquired and run through kanpig on a sequencing experiment (without a `--gqconfig`).
This is the hardest part of GQ calibration. Example VCFs and how they were used are in the
[examples directory](gqcalibration/examples/README.md).

Once the VCF has been genotyped, we reunite the kanpig results with the truth-set genotypes via:

```bash
bcftools merge -m none truth.vcf.gz kanpig.vcf.gz -O z -o merged.vcf.gz
```
Finally, we generate the gqconfig with the command:

```bash
python estimate_params.py merged.vcf.gz my_config
```

This config can now be passed into `kanpig gt --gqconfig my_config.json`

Parameter Details
=================
An important factor for training model parameters is the set of genotypes. Sites with unusually low/high coverage can be
excluded with `--mindp` and `--maxdp`. Furthermore, genotypes from questionable regions of the reference genome should
be excluded by providing a `--bed` file. Generally, the GIAB v1.1 structural variant benchmarking bed file is a fine
default.

As you iterate on different `estimate_params.py` arguments, the most basic way to evaluate if one model is better than
another is by looking at the "Log-likelihood" line in the logging output or the resulting `my_config.json`,
where a higher (less negative) log-likelihood is better. However, you'll also want to look at the QC plots, described
below. For example, the default arguments will only train on already correctly predicted genotypes, which may
over fit the model's parameters. Therefore, it is usually useful to run with `estimate_params.py --all` to also add in the
natural noise of incorrectly genotyped variants, which may give a worse log-likelihood, but may have operate better.

To help avoid over fitting, the `--leaveout` parameter will leave out some fraction of the genotypes (across all states)
from training and then create separate output files for testing.

Output Details
==============

### `gqconfig.json`
This is the main result which can be sent to `kanpig gt --gqconfig` to produce calibrated GTs/GQs.

### `genotypes.csv`
This file holds the information used by `estimate_params.py` with columns of:
* state - Boolean state of if the genotype was predicted correctly
* Ogt - Original truth-set genotype (REF, HET, HOM)
* Kgt - Kanpig predicted genotype (REF, HET, HOM)
* DP - Site's sequencing depth
* AD_ref - Number of reads supporting the reference allele
* AD_alt - Number of reads supporting the alternate allele
* GQ - The default parameters' genotype quality score
* FT - Kanpig's FORMAT/FT field
* KS - Kanpig's FORMAT/KS field
* nKgt - Kanpig's new predicted genotyp (REF, HET, HOM)
* nState - Kanpig's new state
* nGQ - The configured parameters' genotype quality score

If `--leaveout` was used, a second `leaveout.genotypes.csv` will be output.

### Plots

Three main plots are created that can be used to inspect the results. Each plot will have multiple versions.
* The `original` plot uses the input genotype information before parameter estimation.
* The `calibrated` plot uses the genotype information after parameter estimation.
* The `leaveout` plots will have both `original` and `calibrated` versions derived from the left out sites.

#### `ROC.png`
Receiver Operating Curve of genotypes sorted by their GQ (low-to-high).
A ROC that's further towards the top-left is indicative of a more informative GQ, 

#### `GTGQ.png`
This plots the observed genotype accuracy and expected genotype accuracy according to the GQ.
These lines should be highly overlapping when the GQ is well calibrated and actually reflects the probability a genotype
is wrong.

#### `STATE.png`
This plot has three rows corresponding to the REF, HET, and HOM genotype subsets of the data. Each row has four columns.
* Baseline - The distribution of GQs based on their baseline (i.e. true) genotype and colored by their state 
* Kanpig - The distribution of GQs based on their kanpig predicted genotype and colored by their state 
* True GT - The distribution of allele balances (percent of reads supporting the alternate) for correctly predicted
  genotypes
* False GT - The distribution of allele balances (percent of reads supporting the alternate) for incorrectly predicted
  genotypes

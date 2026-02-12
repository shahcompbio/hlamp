# Pipeline to run ECADeMix

This directory contains the Snakemake pipeline that was used to run ECADeMix for the HLAMPs paper, along with required initialization files and parameters. It is primarily intended to be useful for reproducibility purposes, although it can also be used as a template for running ECADeMix on new data.

# Directories/Files

* Snakefile and `rules/` include Snakemake rules referencing scripts in `scripts/`
* `scripts/` includes Python scripts which call functions from the `ecademix` package
* `metadata/` includes the parameters and regions of interest used to generate results for the publication

# Execution

Basic command:
`snakemake run_all_ecademix --configfile config.yaml`

This command can be augmented e.g. with an executor to run in a cluster environment or on the cloud.

# ECADeMix

Python package which uses an MILP to infer amplicons and their multiplicities (loadings) from high-resolution (e.g., 10kb bins) single-cell read counts. The supplied script `infer_amplicons.py` can be used to run ECADeMix. See `pipelines/ECADeMix` for examples scripts which prepare inputs and call functions from this package.

# Setup

`pip install .` in your environment of choice (use of an environment manager like `micromamba` is recommended)

Installation should take no more than a few minutes.

# Dependencies

* `python3`
* `gurobipy` and Gurobi optimizer license
* `numpy`
* `matplotlib`
* `seaborn`
* `scgenome` (https://github.com/mondrian-scwgs/scgenome)
* `pandas`
* `anndata`
* `scipy`
* `scikit-learn`

## Versions used for paper

* `python==3.9.19`
* `gurobipy==12.0.2`
* `numpy==1.23.0`
* `matplotlib==3.9.4`
* `seaborn==0.13.2`
* `scgenome==0.0.19`
* `pandas==2.1.4`
* `anndata==0.10.3`
* `scipy==1.13.1`
* `scikit-learn==1.3.2`

# Input

The main inference function in ECADeMix is `ecademix.solve.optimize_model` which takes as input an `numpy.ndarray` of read counts (dimensions: cells by bins) as well as a variety of optional parameters. For example usage in context, see `infer_amplicons.py`.

The script `infer_amplicons.py` takes as input 3 files: an anndata containing 10kb-binned read counts for regions of interest, an "HMMCopy reads" table containing 500kb copy number profiles and read counts, and an "HMMCopy metrics" table containing cell metrics. The anndata was created using the pipeline in `pipelines/ECADeMix`, and the HMMCopy tables were generated using Mondrian (https://github.com/mondrian-scwgs/mondrian_nf/).

# Usage

Typical runtime varies greatly with the size and complexity of the dataset, as it depends on iterative optimization via Gurobi. Runtime can vary from minutes to several hours.

```bash
python infer_amplicons.py \
    --patient_id <ID> \
    --ilp_input_adata <path> \
    --hmmcopy_reads <path> [--hmmcopy_reads <path> ...] \
    --hmmcopy_metrics <path> [--hmmcopy_metrics <path> ...] \
    --min_quality <float> \
    --max_reads_nonamplicon_cell <float> \
    --n_amplicons <int> \
    --max_amplicon_cn <int> \
    --max_loading <int> \
    --min_loading <int> \
    --lambda_x <float> \
    --lambda_x_l1 <float> \
    --lambda_y <float> \
    --seed <int> \
    --gurobi_license_file <path> \
    --outdir <path>
```

## Example usage for test dataset

This test dataset consists of read counts from 25 cells from patient Lx33. It should run in a few minutes and produce an anndata that matches that in `test/example_output`. Note that you must specify the path to your Gurobi license file as the argument to `--gurobi_license_file`.
```
python infer_amplicons.py \
    --patient_id Lx33 \
    --ilp_input_adata test/input/Lx33_test_adata.h5ad \
    --hmmcopy_reads test/input/Lx33_test_hmmcopy_reads.csv.gz \
    --hmmcopy_metrics test/input/Lx33_test_hmmcopy_metrics.csv.gz \
    --min_quality 0.2 \
    --max_reads_nonamplicon_cell 500 \
    --n_amplicons 1 \
    --max_amplicon_cn 5 \
    --max_loading 50000 \
    --min_loading 0 \
    --lambda_x 100.0 \
    --lambda_x_l1 1.0 \
    --lambda_y 100.0 \
    --plot_data_vmax 350.0 \
    --plot_loading_vmax 100.0 \
    --seed 0 \
    --gurobi_license_file /path/to/gurobi/lic \
    --outdir test/output
```

## Required Arguments

| Argument | Type | Description |
|---|---|---|
| `--patient_id` | string | Patient or sample identifier. Used to select a pre-defined archetype ordering for visualization. |
| `--ilp_input_adata` | path | H5AD file containing the 10 kb single-cell read counts and amplicon annotations (output of upstream pipeline steps). Must have a `var` column `in_amplicon` and a `layers['reads']` entry. |
| `--hmmcopy_reads` | path (repeatable) | CSV file(s) of 500 kb HMMCopy read counts. Pass multiple files by repeating the flag. |
| `--hmmcopy_metrics` | path (repeatable) | CSV file(s) of 500 kb HMMCopy per-cell quality metrics. Pass multiple files by repeating the flag. |
| `--min_quality` | float | Minimum HMMCopy quality score. Cells below this threshold are excluded. |
| `--max_reads_nonamplicon_cell` | float | Maximum total scaled read count used to identify amplicon-free cells. These cells are used to estimate a normalization profile when `--normalize_data` is set. |
| `--n_amplicons` | int | Number of amplicon components to solve for in the ILP. |
| `--max_amplicon_cn` | int | Maximum integer copy number allowed per amplicon bin. |
| `--max_loading` | int | Maximum per-cell loading (amplicon copy dosage) allowed. |
| `--min_loading` | int | Minimum per-cell loading allowed. |
| `--lambda_x` | float | L2 regularization weight on amplicon profiles (X matrix). |
| `--lambda_x_l1` | float | L1 regularization weight on amplicon profiles (X matrix), encouraging sparsity. |
| `--lambda_y` | float | Regularization weight on per-cell loadings (Y matrix). |
| `--seed` | int | Random seed for NMF initialization. |
| `--gurobi_license_file` | path | Path to the Gurobi license file (`gurobi.lic`). Must be set before the Gurobi solver is imported. |
| `--outdir` | path | Directory where output files are written. |

## Optional Arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `--normalize_data` | flag | `False` | If set, divides per-bin read counts by a normalization profile estimated from amplicon-free cells (those with total scaled reads ≤ `--max_reads_nonamplicon_cell`). Also produces a normalization diagnostic plot. |
| `--round_cn_normalize` | flag | `False` | When `--normalize_data` is set, rounds the estimated normalization profile to integer copy numbers before dividing. |
| `--noamplicon_profile_file` | path | `None` | Path to a comma-delimited file containing an externally supplied normalization profile. Applied in addition to the data-derived profile when `--normalize_data` is set. |
| `--manual_init_file` | path | `None` | Path to a comma-delimited file containing a manually specified starting amplicon profile matrix. Overrides NMF initialization. |
| `--cell_weights` | flag | `False` | If set, weights cells during ILP optimization. |
| `--refit_loadings_lambda_y` | float | `NaN` (disabled) | If provided, re-optimizes per-cell loadings after the initial solve using this lambda value, holding amplicon profiles fixed. Useful for capturing more cells with non-zero loadings. |
| `--plot_data_vmax` | float | `500` | Upper limit of the color scale in the residuals data heatmap. |
| `--plot_loading_vmax` | float | `100` | Upper limit of the color scale in the loadings heatmap. |

## Outputs

All outputs are written to `--outdir`.

| File | Description |
|---|---|
| `<patient_id>_result_adata.h5ad` | AnnData containing the ILP solution. `obsm['loadings']` holds per-cell amplicon dosages; `varm['amplicons']` holds inferred per-bin copy-number profiles; `obs['archetype']` holds a patient-specific ordering value for visualization. |
| `<patient_id>_normalization.png` | Diagnostic plot of the normalization profile (only produced when `--normalize_data` is set). |
| Residual / loading heatmaps | QC plots produced by `unmixicon.postprocess.plot_residuals`. |





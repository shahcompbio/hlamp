# eicicle: HLAMP Classifier 


We describe Eicicle, an algorithm to classify HLAMPs as either extrachromosomal (ecDNA) or intrachromosomal (ICamp). Eicicle is based on the following observation. Empirically, the copy number follows a broad distribution when it is driven by the asymmetrical segregation of ecDNA, but a peaked distribution when it is driven by the symmetric segregation of ICamps. Eicicle classifies the cells into those that likely came from the broad ecDNA distribution or the narrow ICamp distribution.


Eicicle takes sets of cells, one per patient, and outputs a classification of each set as coming from either ecDNA or ICamp. The algorithm proceeds as follows.

(1) It first fits a Gaussian mixture to each set of cells, estimating
their distribution.

(2) It then computes two scores from the estimated distributions---an ecDNA score and a mass-in-window score. These scores are indicative of ecDNA and ICamp, respectively.

(3) It labels the sets of cells at the extreme values of these scores
as coming from ecDNA and ICamp.

(4) It runs a linear discriminant analysis to classify the rest of the
sets of cells.



## System requirements

### Software dependencies

Eicicle is a Python package. The following packages are required:

| Package        | Minimum version |
|----------------|-----------------|
| `python`       | 3.10            |
| `numpy`        | 1.24            |
| `pandas`       | 2.0             |
| `scipy`        | 1.10            |
| `scikit-learn` | 1.3             |
| `torch`        | 2.0             |
| `pyro-ppl`     | 1.8             |
| `matplotlib`   | 3.7             |
| `seaborn`      | 0.12            |
| `tqdm`         | 4.65            |
| `pyyaml`       | 6.0             |
| `adjustText`   | 1.0             |
| `jupyter`      | 1.0             |


## Installation guide

### Using Conda

The instructions below install the environment to a custom prefix path. This is useful when the home directory has limited space (as on shared HPC systems). Users with sufficient space in their home directory can replace `--prefix <path>` with `--name eicicle`.

```bash
git clone https://github.com/shahcompbio/hlamp.git
cd hlamp/src/eicicle

# Create the environment at a custom location.
conda env create \
    --prefix /path/to/conda/eicicle \
    --file environment.yml

# Activate it.
conda activate /data1/shahs3/users/salehis/conda/eicicle
```



## Fitting the mixture model

```
python mixture_models_pyro_gaussian.py \
    --dat_path ../../dat/simulated_data.csv.gz \
    --locus CCNE1 \
    --subset_counts -1 \
    --seed 100 \
    --max_M 2 \
    --holdoff_rate 0.01 \
    --outDir . \
    --n_steps 2000 \
    --no_loop \
    --learning_rate .02
```


## Instructions for use

### Input data format

`load_data` (in `mixture_models_pyro_utils.py`) expects a CSV or TSV file (gzip-compressed is supported via the `.csv.gz` / `.tsv.gz` extension) with:

- A `cell_id` column identifying each cell.
- One column per locus, named by the locus identifier passed via `--locus`. Values are copy numbers (non-negative; NaNs are dropped).

Example:

```
cell_id,ARTLocus
cell_0,15
cell_1,6
cell_2,4
...
```

By default, `load_data` drops NaNs, drops negative values, and (if `do_prune_left=True`) prunes spurious low-copy-number tails. See the function docstring for the full set of options.


## Multi Sample Mode

`ecDNA-score` can be computed in the multi-case mode in which both positive and negative cases are observed. 
In the multi-case mode we compute the Z-score for mean and s.d. of the components with the largest mixing proportion, then rescale the resulting quantity to lie between 0 and 1. 
We report this mode in our manuscript. 



## Single Sample Mode

It is possible to compute the `ecDNA-score` in a single-sample mode. 
In this mode, we do not Z-score the mean and s.d. of the major component, nor do we rescale it to lie between 0 and 1. 



## Files

```
mixture_models_pyro_utils.py    Utility functions for loading and filtering data.
mixture_models_pyro_gaussian.py Implementation of the Gaussian mixture model (entry point).
mixture_models_pyro_model.py    Pyro model definitions used by the fitter.
mixture_models_pyro_fits.py     Helpers for exporting fit results.
mixture_models_pyro_globals.py  Global constants and width helpers.
mixture_utils.py                Helpers for writing YAML and picking the best model.
compute_bayes_scores_utils.py   Utility functions for computing the summary scores.
compute_bayes_scores.ipynb      Notebook that computes the summary scores from model inference results.
compute_bayes_scores_single.ipynb  Single-sample variant of the above.
decision_boundary.ipynb         Uses LDA to fit a classifier on the data.
environment.yml                 Conda environment specification.
requirements.txt                Pip requirements file (alternative to conda).
```



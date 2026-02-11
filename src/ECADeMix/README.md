# ECADeMix

Python package which uses an MILP to infer amplicons and their multiplicities (loadings) from high-resolution (e.g., 10kb bins) single-cell read counts. See `pipelines/ECADeMix` for example scripts which call functions from this package.

# Setup

`pip install .`

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

* `python=3.9.19`
* `gurobipy==12.0.2`
* `numpy==1.23.0`
* `matplotlib==3.9.4`
* `seaborn==0.13.2`
* `scgenome==0.0.19`
* `pandas==2.1.4`
* `anndata==0.10.3`
* `scipy==1.13.1`
* `scikit-learn==1.3.2`



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



## Files

```
mixture_models_pyro_utils.py: Utility functions for loading and filtering data.mixture_models_pyro_gaussian.py: Implementation of the Gaussian mixture model.
compute_bayes_scores.ipynb: Notebook that computes the summary scores from model inference results.
compute_bayes_scores_utils.py: Utility functions for computing the summary scores.
decision_boundary.ipynb: Uses LDA to fit a classifier on the data. 
```


## Fitting the mixture model

```
python mixture_models_pyro_gaussian.py \
    --dat_path path_to_data \
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


## Multi Sample Mode

`ecDNA-score` can be computed in the multi-case mode in which both positive and negative cases are observed. 
In the multi-case mode we compute the Z-score for mean and s.d. of the components with the largest mixing proportion, then rescale the resulting quantity to lie between 0 and 1. 
We report this mode in our manuscript. 



## Single Sample Mode

It is possible to compute the `ecDNA-score` in a single-sample mode. 
In this mode, we do not Z-score the mean and s.d. of the major component, nor do we rescale it to lie between 0 and 1. 




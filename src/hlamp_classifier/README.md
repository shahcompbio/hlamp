## Files

```
mixture_models_pyro_utils.py: Utility functions for loading and filtering data.mixture_models_pyro_gaussian.py: Implementation of the Gaussian mixture model.
compute_bayes_scores.ipynb: Notebook that computes the summary scores from model inference results.
compute_bayes_scores_utils.py: Utility functions for computing the summary scores.
decision_boundary.ipynb: Uses LDA to fit a classifier on the data. 
```


## Fitting the mixture model

```
python src/mixture_model/mixture_models_pyro_gaussian.py \
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
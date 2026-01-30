## Files

```
bfb_simul.py: Main simulation code with CLI. bfb_simul_viz.py: Utility functions for visualization in form of histograms.
```


## Example run

```
python src/bfb_simul.py \
    --seed 20001 \
    --N_CELLS 100000000 \
    --N_GENERATION 5000 \
    --N_INITIAL_CELLS 1 \
    --n_nDNA 2 \
    --death_rate 0.01 \
    --div_rate_base 1.0 \
    --p_bfb 1.0 \
    --p_damaged 0.0 \
    --p_fix 0.00 \
    --p_wgd 0.00 \
    --damaged_death_coeff 1 \
    --division_mode logistic \
    --damaged_div_mult 0.0 \
    --damage_bfb_multiplier 1 \
    --p_ecDNA_segregation 0.6 \
    --death_scenario low_and_high \
    --death_rate_params 'death_low_mult=2.0,death_mid_mult=0.001,death_high_mult=4.0' \
    --outDir ./results/bfb_simul/
    ```
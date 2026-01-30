#!/bin/python
# Run the stochastic agent based simulator for copy number change.

import os
import pickle
import yaml
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import argparse
from bfb_simul_viz import plot_population_size, plot_histogram, plot_subsample_histogram


def parse_death_rate_overrides(death_rate_params):
    overrides = {}
    for item in death_rate_params.split(","):
        if "=" in item:
            key, value = item.split("=")
            try:
                value = float(value)
                if value.is_integer():
                    value = int(value)
            except ValueError:
                pass  # keep as string if not a number
            overrides[key.strip()] = value
    return overrides


def get_death_rate_scenario(scenario_name, **overrides):
    """
    Return kwargs for get_death_rate() based on a named scenario.

    Override any returned parameter by passing it as **overrides.
    """

    # Define fixed values for the low and high cutoffs
    DEFAULT_LOW_CUT = 100
    DEFAULT_HIGH_CUT = 500
    # DEFAULT_LOW_CUT = 10
    # DEFAULT_HIGH_CUT = 100

    scenario_name = str(scenario_name).strip().lower()

    if scenario_name in {"neutral", "no_penalty", "flat"}:
        params = {
            "mode": "neutral",
        }

    elif scenario_name in {"low_and_high", "both_sides", "u_shaped"}:
        params = {
            "mode": "cn_thresholds",
            "low_cut": DEFAULT_LOW_CUT,
            "high_cut": DEFAULT_HIGH_CUT,
            "death_low_mult": 2.0,
            # "death_mid_mult": 10.0,
            # "death_mid_mult": .1,
            "death_mid_mult": 0.001,
            # "death_high_mult": 3.0,
            "death_high_mult": 4.0,
        }

    elif scenario_name in {"high_only", "penalize_high"}:
        params = {
            "mode": "cn_thresholds",
            "low_cut": DEFAULT_LOW_CUT,
            "high_cut": DEFAULT_HIGH_CUT,
            "death_low_mult": 1.0,  # no low penalty
            "death_mid_mult": 1.0,
            "death_high_mult": 3.0,  # high penalty only
        }
    # one for mid_and_high

    elif scenario_name in {"mid_and_high", "medium_and_high"}:
        params = {
            "mode": "cn_thresholds",
            "low_cut": DEFAULT_LOW_CUT,
            "high_cut": DEFAULT_HIGH_CUT,
            "death_low_mult": 1.0,  # no low penalty
            "death_mid_mult": 3.0,  # medium penalty
            "death_high_mult": 3.0,  # high penalty
        }

    elif scenario_name in {"carrying_capacity", "capacity", "density"}:
        params = {
            "mode": "carrying_capacity",
            # death is multiplied by (1 + cc_strength * (n_cells/N_CELLS)**cc_power)
            "cc_strength": 5.0,
            "cc_power": 1.0,
        }

    else:
        raise ValueError(
            f"Unknown death-rate scenario: {scenario_name}. "
            f"Valid: neutral, low_and_high, high_only, carrying_capacity"
        )
    params.update(overrides)
    return params


def get_default_params():
    # ---------------------------
    # Constants
    # ---------------------------
    param_dict = {}
    param_dict["seed"] = 100001  # Random seed
    param_dict["N_CELLS"] = int(1e7)  # Maximum allowed cells
    param_dict["N_GENERATION"] = 1000  # Maximum generations
    param_dict["N_INITIAL_CELLS"] = 10  # Start with one cell
    param_dict["n_nDNA"] = 2  # Initial copies per cell
    param_dict["n_ecDNA"] = 0  # Initial copies per cell
    param_dict["death_rate"] = 0.001  # Base death rate
    param_dict["div_rate_base"] = 0.9  # Base division rate
    param_dict["p_reintegrate"] = 0.000  # ecDNA reintegration probability
    param_dict["p_n_eDNA"] = 0.0  # ecDNA creation probability (from nDNA)
    param_dict["p_bfb"] = 0.5  # BFB event probability
    param_dict["p_damaged"] = 0.9  # Damage probability given a BFB event
    param_dict["p_fix"] = 0.001  # Probability that damaged DNA is fixed
    param_dict["p_wgd"] = 0.0001  # Whole genome duplication probability
    param_dict["damaged_death_coeff"] = 50  # Death rate multiplier for damaged cells
    param_dict["p_ecDNA_segregation"] = 0.6
    return param_dict


def get_division_rate(
    ndna, N_CELLS, div_rate_base=0.9, division_mode="logistic", **kwargs
):
    """
    Compute cell-specific division probabilities based on the overall population size
    and the cell's nuclear DNA (ndna) copy number, applying different multipliers based
    on two thresholds, N1 and N2.

    Parameters
    ----------
    ndna : np.array
        Array of nuclear DNA copy numbers for each cell.
    N_CELLS : int
        Maximum allowed number of cells.
    div_rate_base : float
        Base division rate.
    division_mode : str
        Division mode ('logistic', 'linear', 'exponential', or 'constant').
    **kwargs :
        N1 : int, optional
            Threshold for ndna below which multiplier N1_m is applied (default: 10).
        N2 : int, optional
            Threshold for ndna above which multiplier N1_m2 is applied (default: 20).
        N1_m : float, optional
            Multiplier for cells with ndna < N1 (default: 0.8).
        N1_m1 : float, optional
            Multiplier for cells with N1 <= ndna <= N2 (default: 1.0).
        N1_m2 : float, optional
            Multiplier for cells with ndna > N2 (default: 1.2).

    Returns
    -------
    div_prob : np.array
        Array of division probabilities for each cell.
    """
    # Global division probability based on overall cell count.
    n_cells = ndna.shape[0]
    if division_mode == "logistic":
        sharpness_factor = 1  # Adjust this value if needed
        global_div = div_rate_base * (
            1 - sharpness_factor * (np.log(n_cells + 1) / np.log(N_CELLS + 1))
        )
    elif division_mode == "linear":
        global_div = div_rate_base * (1 - n_cells / N_CELLS)
    elif division_mode == "exponential":
        global_div = div_rate_base * np.exp(-n_cells / N_CELLS)
    elif division_mode == "constant":
        global_div = div_rate_base
    else:
        raise ValueError(f"Unsupported division mode: {division_mode}")

    # Return an array of the size ndna with the global_div value
    return global_div * np.ones_like(ndna)

    # Get threshold values and multipliers from kwargs.
    low_cut = kwargs.get("N1", 50)
    high_cut = kwargs.get("N2", 100)
    # The multipliers for the division probability.
    low_mult = kwargs.get("N1_m", 0.7)
    mid_mult = kwargs.get("N1_m1", 1.0)
    high_mult = kwargs.get("N1_m2", 0.9)

    # Compute cell-specific multiplier based on ndna thresholds.
    multiplier = np.where(
        ndna < low_cut, low_mult, np.where(ndna > high_cut, high_mult, mid_mult)
    )

    # Final division probability for each cell.
    div_prob = global_div * multiplier
    div_prob = np.clip(div_prob, 0, 1)

    return div_prob


# def get_death_rate(ndna, base_death_rate, **kwargs):
#     """
#     Compute cell-specific death probabilities based on the cell's nuclear DNA (ndna) copy number,
#     applying different multipliers according to specified thresholds. Cells with lower nDNA
#     receive a higher multiplier, thereby increasing their death probability.

#     Parameters
#     ----------
#     ndna : np.array
#         Array of nuclear DNA copy numbers for each cell.
#     base_death_rate : float
#         Base death rate.
#     **kwargs :
#         low_cut : int, optional
#             Threshold for nDNA below which a higher death multiplier is applied (default: 10).
#         high_cut : int, optional
#             Threshold for nDNA above which a lower death multiplier is applied (default: 50).
#         death_low_mult : float, optional
#             Multiplier for cells with nDNA < low_cut (default: 1.2).
#         death_mid_mult : float, optional
#             Multiplier for cells with low_cut <= nDNA <= high_cut (default: 1.0).
#         death_high_mult : float, optional
#             Multiplier for cells with nDNA > high_cut (default: 0.8).

#     Returns
#     -------
#     death_rate_arr : np.array
#         Array of cell-specific death rates after applying the nDNA-based multipliers.
#     """
#     #return base_death_rate
#     low_cut = kwargs.get("low_cut", 50)
#     high_cut = kwargs.get("high_cut", 100)
#     death_low_mult = kwargs.get("death_low_mult", 1.5)
#     death_mid_mult = kwargs.get("death_mid_mult", 1.0)
#     death_high_mult = kwargs.get("death_high_mult", 1.1)

#     # Determine the multiplier for each cell based on its nDNA:
#     multiplier = np.where(
#         ndna < low_cut,
#         death_low_mult,
#         np.where(ndna > high_cut, death_high_mult, death_mid_mult),
#     )

#     death_rate_arr = base_death_rate * multiplier
#     death_rate_arr = np.clip(death_rate_arr, 0, 1)

#     return death_rate_arr


def get_death_rate(ndna, base_death_rate, N_CELLS=None, mode="cn_thresholds", **kwargs):
    """
    Compute cell-specific death probabilities.

    Supported modes
    ---------------
    mode="neutral":
        death_rate_arr = base_death_rate (no CN penalty)
    mode="cn_thresholds":
        Apply low/mid/high multipliers based on low_cut/high_cut
    mode="carrying_capacity":
        Multiply base_death_rate by a density factor that depends on n_cells / N_CELLS
        Optionally, you can still pass cn_threshold params too, but by default it's pure density.
    """
    n_cells = ndna.shape[0]

    if mode == "neutral":
        death_rate_arr = base_death_rate * np.ones_like(ndna, dtype=float)

    elif mode == "cn_thresholds":
        low_cut = kwargs.get("low_cut", 50)
        high_cut = kwargs.get("high_cut", 100)
        death_low_mult = kwargs.get("death_low_mult", 1.5)
        death_mid_mult = kwargs.get("death_mid_mult", 1.0)
        death_high_mult = kwargs.get("death_high_mult", 1.1)

        multiplier = np.where(
            ndna < low_cut,
            death_low_mult,
            np.where(ndna > high_cut, death_high_mult, death_mid_mult),
        )
        death_rate_arr = base_death_rate * multiplier

    elif mode == "carrying_capacity":
        if N_CELLS is None:
            raise ValueError("mode='carrying_capacity' requires N_CELLS to be provided")

        cc_strength = kwargs.get("cc_strength", 5.0)
        cc_power = kwargs.get("cc_power", 1.0)

        density = n_cells / float(N_CELLS)
        density_factor = 1.0 + cc_strength * (density**cc_power)

        death_rate_arr = (
            base_death_rate * density_factor * np.ones_like(ndna, dtype=float)
        )

    else:
        raise ValueError(f"Unsupported death-rate mode: {mode}")

    return np.clip(death_rate_arr, 0, 1)


def vectorized_narrow_split(totals, mean_frac=0.5, sd_frac=0.05):
    """
    Splits an array of totals into two parts using a low-variance random approach.

    Parameters
    ----------
    totals : np.array
        Array of total values to be split.
    mean_frac : float, optional
        Fraction of the total to use as the mean for the split (default: 0.5).
    sd_frac : float, optional
        Fraction of the total to use as the standard deviation for the split (default: 0.05).

    Returns
    -------
    split_vals : np.array
        Array of split values (for daughter1) obtained by sampling from a narrow normal distribution,
        then clipping to [0, total] and rounding to the nearest integer.
    """
    means = totals * mean_frac
    sds = totals * sd_frac
    # Vectorized normal samples; note that some samples may fall outside [0, total].
    splits = np.random.normal(loc=means, scale=sds)
    # Clip to ensure the split is between 0 and the total
    splits = np.clip(splits, 0, totals)
    # Round to the nearest integer (and cast to integer type)
    return np.rint(splits).astype(np.int64)


def simulate(
    seed,
    N_CELLS,
    N_GENERATION,
    N_INITIAL_CELLS,
    n_nDNA,
    n_ecDNA,
    death_rate,
    div_rate_base,
    p_reintegrate,
    p_n_eDNA,
    p_bfb,
    p_damaged,
    p_fix,
    p_wgd,
    damaged_death_coeff,
    division_mode="logistic",
    damaged_div_mult=0.01,
    damage_bfb_multiplier=2.0,
    call_back_func=None,
    percent_complete=0.999,
    p_ecDNA_segregation=0.6,
    death_scenario="low_and_high",
    death_rate_params=None,
    symmetric_division=False,
):
    """
    Run the stochastic agent based simulator for copy number change.

    Parameters
    ----------
    seed : int
        Random seed.
    N_CELLS : int
        Maximum allowed number of cells.
    N_GENERATION : int
        Maximum number of generations.
    N_INITIAL_CELLS : int
        Number of identical cells at the start.
    n_nDNA : int
        Initial copies of nuclear DNA (nDNA) per cell.
    n_ecDNA : int
        Initial copies of extrachromosomal DNA (ecDNA) per cell.
    death_rate : float
        Base death rate.
    div_rate_base : float
        Base division rate.
    p_reintegrate : float
        ecDNA reintegration probability.
    p_n_eDNA : float
        ecDNA creation probability (from nDNA).
    p_bfb : float
        Base BFB event probability.
    p_damaged : float
        Damage probability given a BFB event.
    p_fix : float
        Probability that damaged DNA is fixed.
    p_wgd : float
        Whole genome duplication probability.
    damaged_death_coeff : float
        Death rate multiplier for damaged cells.
    division_mode : str, optional
        Division mode ('logistic', 'linear', 'exponential', or 'constant').
    call_back_func : function, optional
        Called at the end of each generation.
    p_ecDNA_segregation : float
        Probability of ecDNA segregation to daughter cell. (.5 equal, >0.5 biased)
    death_rate_params : string, optional
        Death rate parameters in the format key1=value1;key2=value2;...
    symmetric_division : bool
        If True, all divisions are symmetric (equal split of nDNA and ecDNA).

    Returns
    -------
    summary : np.array
        Population summary as (generation, number of cells, total nDNA, total ecDNA).
    final_ndna : np.array
        Final nDNA copy numbers for all cells.
    final_ecdna : np.array
        Final ecDNA copy numbers for all cells.
    final_damaged : np.array
        Final damage status for all cells.
    """
    assert (
        p_ecDNA_segregation >= 0 and p_ecDNA_segregation <= 1
    ), "p_ecDNA_segregation must be between 0 and 1"
    np.random.seed(seed)
    ndna = np.full(N_INITIAL_CELLS, n_nDNA, dtype=np.int64)
    ecdna = np.full(N_INITIAL_CELLS, n_ecDNA, dtype=np.int64)
    damaged = np.zeros(N_INITIAL_CELLS, dtype=bool)

    # Parse death_rate_params if provided
    overrides = (
        {}
        if death_rate_params is None
        else parse_death_rate_overrides(death_rate_params)
    )

    summary = []
    pbar = tqdm(enumerate(range(N_GENERATION)))
    n_cells_at_iter = []
    STOP_STATUS = ""
    for ii, gen in pbar:
        n0 = ndna.shape[0]  # Number of cells at the start of the generation.
        if n0 >= N_CELLS:
            print(f"\nReached maximum number of cells at generation {gen}")
            STOP_STATUS = "max_cells_reached"
            break
        if n0 == 0:
            print(f"\nAll cells are dead at generation {gen}")
            STOP_STATUS = "all_cells_dead"
            break

        # ---------------------------
        # Damage Fixing
        # ---------------------------
        fix_rand = np.random.rand(n0)
        damaged = np.where(damaged & (fix_rand < p_fix), False, damaged)

        # ---------------------------
        # Death Phase
        # ---------------------------
        # base_death_rates = get_death_rate(
        #     ndna,
        #     death_rate,
        #     # low_cut=10,
        #     # high_cut=50,
        #     low_cut=100,
        #     high_cut=500,
        #     death_low_mult=2.0,
        #     death_mid_mult=.8,
        #     #death_high_mult=0.8,
        #     death_high_mult=3.0,
        # )
        death_params = get_death_rate_scenario(
            scenario_name=death_scenario, **(overrides)
        )

        base_death_rates = get_death_rate(
            ndna,
            death_rate,
            N_CELLS=N_CELLS,
            **death_params,
        )

        # If cells are damaged, apply the damaged_death_coeff.
        current_death_rate = np.where(
            damaged, base_death_rates * damaged_death_coeff, base_death_rates
        )
        current_death_rate = np.clip(current_death_rate, 0, 1)
        death_rand = np.random.rand(n0)
        survive_mask = (death_rand >= current_death_rate) & (ndna >= 2)

        ndna = ndna[survive_mask]
        ecdna = ecdna[survive_mask]
        damaged = damaged[survive_mask]

        # ---------------------------
        # Division Decision
        # ---------------------------
        div_prob = get_division_rate(
            ndna,
            N_CELLS,
            div_rate_base=div_rate_base,
            division_mode=division_mode,
            n_nDNA=n_nDNA,
        )
        # Make sure that damaged cells divide more slowly.
        # TODO: set the damage div rate
        # div_prob[damaged] = 0.001
        div_prob[damaged] = div_prob[damaged] * damaged_div_mult
        div_rand = np.random.rand(ndna.shape[0])
        dividing_mask = div_rand <= div_prob
        nondividing_mask = ~dividing_mask

        ndna_nondiv = ndna[nondividing_mask]
        ecdna_nondiv = ecdna[nondividing_mask]
        damaged_nondiv = damaged[nondividing_mask]

        # ---------------------------
        # Division (for dividing cells)
        # ---------------------------
        ndna_div = ndna[dividing_mask]
        ecdna_div = ecdna[dividing_mask]
        damaged_div = damaged[dividing_mask]
        n_div = ndna_div.shape[0]

        if n_div > 0:
            # --- Duplication ---
            wgd_rand = np.random.rand(n_div)
            new_ndna = np.where(wgd_rand < p_wgd, ndna_div * 4, ndna_div * 2)
            new_ecdna = np.where(wgd_rand < p_wgd, ecdna_div * 4, ecdna_div * 2)

            # --- BFB and DNA Damage ---
            # Increase the BFB probability for cells that are already damaged.
            effective_p_bfb = np.where(
                damaged_div, p_bfb * damage_bfb_multiplier, p_bfb
            )  # Damaged cells are more likely to undergo BFB.
            effective_p_bfb = np.clip(effective_p_bfb, 0, 1)

            bfb_rand = np.random.rand(n_div)
            is_bfb = bfb_rand < effective_p_bfb

            damage_rand = np.random.rand(n_div)
            new_damaged = np.where(
                is_bfb & (damage_rand < p_damaged), True, damaged_div
            )

            # --- nDNA Splitting ---
            # ndna_div_daughter1 = np.empty(n_div, dtype=np.int64)
            # idx_bfb = np.where(is_bfb)[0]
            # if idx_bfb.size > 0:
            #     ndna_div_daughter1[idx_bfb] = np.random.binomial(
            #         new_ndna[idx_bfb], p_ecDNA_segregation
            #     )
            # idx_nonbfb = np.where(~is_bfb)[0]
            # if idx_nonbfb.size > 0:
            #     ndna_div_daughter1[idx_nonbfb] = new_ndna[idx_nonbfb] // 2
            # ndna_div_daughter2 = new_ndna - ndna_div_daughter1
            ####################################################
            ndna_div_daughter1 = np.empty(n_div, dtype=np.int64)
            if symmetric_division:
                ndna_div_daughter1[:] = new_ndna // 2
            else:
                idx_bfb = np.where(is_bfb)[0]
                if idx_bfb.size > 0:
                    ndna_div_daughter1[idx_bfb] = np.random.binomial(
                        new_ndna[idx_bfb], p_ecDNA_segregation
                    )
                idx_nonbfb = np.where(~is_bfb)[0]
                if idx_nonbfb.size > 0:
                    ndna_div_daughter1[idx_nonbfb] = new_ndna[idx_nonbfb] // 2
            ndna_div_daughter2 = new_ndna - ndna_div_daughter1
            ####################################################

            # Set daughter cell states.
            d1_ndna = ndna_div_daughter1
            d2_ndna = ndna_div_daughter2
            M_d1 = np.zeros(n_div, dtype=np.int64)
            M_d2 = np.zeros(n_div, dtype=np.int64)
            d1_ecdna = M_d1
            d2_ecdna = M_d2
            d1_damaged = new_damaged
            d2_damaged = new_damaged

            ndna_divided = np.concatenate([d1_ndna, d2_ndna])
            ecdna_divided = np.concatenate([d1_ecdna, d2_ecdna])
            damaged_divided = np.concatenate([d1_damaged, d2_damaged])
        else:
            ndna_divided = np.empty(0, dtype=np.int64)
            ecdna_divided = np.empty(0, dtype=np.int64)
            damaged_divided = np.empty(0, dtype=bool)

        # ---------------------------
        # New Generation
        # ---------------------------
        ndna = np.concatenate([ndna_nondiv, ndna_divided])
        ecdna = np.concatenate([ecdna_nondiv, ecdna_divided])
        damaged = np.concatenate([damaged_nondiv, damaged_divided])
        if ndna.size == 0:
            print(f"All cells are dead at generation {gen}")
            STOP_STATUS = "all_cells_dead"
            break
        total_ndna = ndna.sum()
        total_ecdna = ecdna.sum()
        summary.append((gen + 1, ndna.size, total_ndna, total_ecdna))
        pct_complete = max((gen + 1) / N_GENERATION, ndna.size / N_CELLS)
        max_n_dna = np.max(ndna)
        if pct_complete > percent_complete:
            print(f"\nEarly stopping at generation {gen}")
            STOP_STATUS = "early_stopping"
            break
        pbar.set_description(
            f"Generation: {gen+1}/{N_GENERATION}, Cells: {ndna.size}/{N_CELLS}, Pct: {pct_complete:.2f}, MaxCN: {max_n_dna}, nDNA: {total_ndna}, ecDNA: {total_ecdna}"
        )
        if call_back_func is not None:
            call_back_func(ndna, ecdna, damaged, gen, stop_status_str=STOP_STATUS)

        n_cells_at_iter.append(ndna.size)

        # Check if that in the last 50 generations the number of cells has not changed much
        # if len(n_cells_at_iter) > 50:
        #     if np.std(n_cells_at_iter[-50:]) < 100:
        #         print(f"Early stopping at generation {gen}")
        #         break
        # every 10 generations, print the std
        # if gen % 10 == 0:
        #     print(f"std of the last 50 generations: {np.std(n_cells_at_iter[-50:])}")

    return np.array(summary), ndna, ecdna, damaged


def get_rnd_dirname():
    """
    Generate a random filename. Adds in data and time
    """
    import string
    import random
    import datetime

    now = datetime.datetime.now()
    date_time = now.strftime("%Y-%m-%d-%H-%M-%S")
    return (
        "".join(random.choices(string.ascii_lowercase + string.digits, k=8)) + date_time
    )


def main(outDir=None, loadResults=False, **kwargs):
    """
    Run the stochastic agent based simulator for copy number change.

    Runs the simulation and plots:
    - Number of cells over generations
    - Number of nDNA copies over generations
    - Number of ecDNA copies over generations
    - Number of damaged cells over generations
    """
    # Setup the output directory and plotting filenames
    if outDir is not None:
        os.makedirs(outDir, exist_ok=True)
        # Add rnd subdir for this run
        outDir = os.path.join(outDir, get_rnd_dirname())
        # Configurations path
        conf_path = os.path.join(outDir, "config.yaml")
        # Results paths
        resPath1 = os.path.join(outDir, "summary.pkl")
        resTmpDir = os.path.join(outDir, "iters")
        os.makedirs(resTmpDir, exist_ok=True)
        # Figure paths
        figuresDir = os.path.join(outDir, "figures")
        os.makedirs(figuresDir, exist_ok=True)
        figPath1 = os.path.join(figuresDir, "population_size.png")
        figPath2 = os.path.join(figuresDir, "histogram.png")
        figPath3 = os.path.join(figuresDir, "subsample_histogram.png")
    else:
        figPath1 = None
        figPath2 = None
        figPath3 = None
    if loadResults:
        with open(resPath1, "rb") as f:
            results_dict = pickle.load(f)
        summary = results_dict["summary"]
        final_ndna = results_dict["final_ndna"]
        final_ecdna = results_dict["final_ecdna"]
        final_damaged = results_dict["final_damaged"]
        plot_population_size(summary=summary, final_ndna=final_ndna, figPath=figPath1)
        plot_histogram(final_ndna=final_ndna, figPath=figPath2)
        plot_subsample_histogram(final_ndna=final_ndna, figPath=figPath3)
        print(f"Results saved to {outDir}")
        return

    # ---------------------------
    # Save the configuration
    # ---------------------------
    with open(conf_path, "w") as f:
        yaml.dump(kwargs, f)

    # ---------------------------
    # Run the Simulation
    # ---------------------------
    def call_back_saver(
        ndna,
        ecdna,
        damaged,
        gen,
        period=100,
        do_plot=True,
        save_pickle=True,
        stop_status_str=None,
    ):
        """
        Save the results every KK generations
        """
        if (
            gen % period == 0
            or gen == kwargs["N_GENERATION"] - 1
            or (stop_status_str is not None and stop_status_str == "early_stopping")
        ):
            results_dict = {
                "final_ndna": ndna,
                "final_ecdna": ecdna,
                "final_damaged": damaged,
                "generation": gen,
                "stop_status_str": stop_status_str,
            }
            if save_pickle:
                resPathTmp = os.path.join(resTmpDir, f"summary_{gen}.pkl")
                with open(resPathTmp, "wb") as f:
                    pickle.dump(results_dict, f)
            if do_plot:
                plot_histogram(
                    final_ndna=ndna,
                    figPath=os.path.join(resTmpDir, f"histogram_{gen}.png"),
                )
                plot_subsample_histogram(
                    final_ndna=ndna,
                    figPath=os.path.join(resTmpDir, f"subsample_histogram_{gen}.png"),
                )

    kwargs["call_back_func"] = call_back_saver
    summary, final_ndna, final_ecdna, final_damaged = simulate(**kwargs)
    # ---------------------------
    # Save the results
    # ---------------------------
    results_dict = {
        "summary": summary,
        "final_ndna": final_ndna,
        "final_ecdna": final_ecdna,
        "final_damaged": final_damaged,
    }
    with open(resPath1, "wb") as f:
        pickle.dump(results_dict, f)
    # ---------------------------
    # Plot the Results
    # ---------------------------
    plot_population_size(summary=summary, final_ndna=final_ndna, figPath=figPath1)
    plot_histogram(final_ndna=final_ndna, figPath=figPath2)
    plot_subsample_histogram(final_ndna=final_ndna, figPath=figPath3)
    print(f"Results saved to {outDir}")


if __name__ == "__main__":
    """
    @Example:

    python bfb_simul.py \
    --seed 100 \
    --N_CELLS 100000 \
    --N_GENERATION 200 \
    --N_INITIAL_CELLS 10 \
    --n_nDNA 2 \
    --n_ecDNA 0 \
    --death_rate 0.001 \
    --div_rate_base 0.9 \
    --p_reintegrate 0.000 \
    --p_n_eDNA 0.0 \
    --p_bfb 0.5 \
    --p_damaged 0.9 \
    --p_fix 0.001 \
    --p_wgd 0.0001 \
    --damaged_death_coeff 50 \
    --division_mode logistic \
    --outDir ../../results/bfb_simul_feb_07

    python bfb_simul.py \
    --seed 20001 \
    --N_CELLS 10000000 \
    --N_GENERATION 5000 \
    --N_INITIAL_CELLS 10 \
    --n_nDNA 2 \
    --n_ecDNA 0 \
    --death_rate 0.001 \
    --div_rate_base 1.0 \
    --p_bfb 0.5 \
    --p_damaged 0.9 \
    --p_fix 0.0001 \
    --p_wgd 0.00 \
    --damaged_death_coeff 10 \
    --division_mode linear \
    --damaged_div_mult 0.01 \
    --damage_bfb_multiplier 3 \
    --outDir ../../results/bfb_simul_feb_07


    """
    parser = argparse.ArgumentParser(
        description="Run the stochastic agent based simulator for copy number change. "
    )
    parser.add_argument("--seed", type=int, default=100001, help="Random seed")
    parser.add_argument(
        "--N_CELLS", type=int, default=int(1e7), help="Maximum allowed cells"
    )
    parser.add_argument(
        "--N_GENERATION", type=int, default=1000, help="Maximum generations"
    )
    parser.add_argument(
        "--N_INITIAL_CELLS", type=int, default=10, help="Start with one cell"
    )
    parser.add_argument("--n_nDNA", type=int, default=2, help="Initial copies per cell")
    parser.add_argument(
        "--n_ecDNA", type=int, default=0, help="Initial copies per cell"
    )
    parser.add_argument(
        "--death_rate", type=float, default=0.001, help="Base death rate"
    )
    parser.add_argument(
        "--div_rate_base", type=float, default=0.9, help="Base division rate"
    )
    parser.add_argument(
        "--p_reintegrate",
        type=float,
        default=0.000,
        help="ecDNA reintegration probability",
    )
    parser.add_argument(
        "--p_n_eDNA",
        type=float,
        default=0.0,
        help="ecDNA creation probability (from nDNA)",
    )
    parser.add_argument(
        "--p_bfb", type=float, default=0.5, help="BFB event probability"
    )
    parser.add_argument(
        "--p_damaged",
        type=float,
        default=0.9,
        help="Damage probability given a BFB event",
    )
    parser.add_argument(
        "--p_fix",
        type=float,
        default=0.001,
        help="Probability that damaged DNA is fixed",
    )
    parser.add_argument(
        "--p_wgd",
        type=float,
        default=0.0001,
        help="Whole genome duplication probability",
    )
    parser.add_argument(
        "--damaged_death_coeff",
        type=float,
        default=50,
        help="Death rate multiplier for damaged cells",
    )
    parser.add_argument(
        "--division_mode",
        type=str,
        default="logistic",
        help="Division mode ('logistic', 'linear', 'exponential', 'constant')",
    )
    parser.add_argument(
        "--outDir",
        type=str,
        default=None,
        help="Output directory for the simulation results",
    )
    # add a boolean for load the results
    parser.add_argument(
        "--loadResults",
        action="store_true",
        help="Load the results from a previous run",
    )
    parser.add_argument(
        "--damaged_div_mult",
        type=float,
        default=0.001,
        help="Division rate multiplier for damaged cells",
    )
    parser.add_argument(
        "--damage_bfb_multiplier",
        type=float,
        default=2.0,
        help="BFB event probability multiplier for damaged cells",
    )
    # add a param for percent complete early stopping
    parser.add_argument(
        "--percent_complete",
        type=float,
        default=0.99,
        help="Percent complete early stopping",
    )
    # add a probability of segregation of ecDNA
    parser.add_argument(
        "--p_ecDNA_segregation",
        type=float,
        default=0.6,
        help="Probability of ecDNA segregation to daughter cell",
    )
    parser.add_argument(
        "--death_scenario",
        type=str,
        default="low_and_high",
        help="Death-rate scenario: neutral, low_and_high, high_only, carrying_capacity",
    )
    # Add a stirng argument that takes the death rates (d1;d2;d3) --> d1, d2, d3
    parser.add_argument(
        "--death_rate_params",
        type=str,
        default="",
        help="Death rate parameters in the format 'key1=value1;key2=value2;...'",
    )
    parser.add_argument(
        "--symmetric_division",
        action="store_true",
        help="Force symmetric splitting on division (override binomial BFB split)",
    )
    args = parser.parse_args()
    main(**vars(args))

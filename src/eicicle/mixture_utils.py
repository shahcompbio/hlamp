import yaml
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def convert_numpy_to_python(obj):
    """Recursively convert any numpy arrays to lists."""
    if isinstance(obj, dict):
        return {k: convert_numpy_to_python(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_to_python(x) for x in obj]
    elif isinstance(obj, tuple):
        return tuple(convert_numpy_to_python(x) for x in obj)
    elif isinstance(obj, np.ndarray):
        # Convert arrays to (nested) Python lists
        return obj.tolist()
    elif isinstance(obj, np.generic):
        # Convert NumPy scalars to native Python scalars (float, int)
        return obj.item()
    else:
        return obj


def write_results_to_yaml(results, outPath, verbose=False):
    # Convert torch tensors to numpy arrays
    results = convert_numpy_to_python(results)
    if verbose:
        print(results)
    # Don't use numpy binary, just plain text
    with open(outPath, "w") as f:
        yaml.dump(
            results,
            f,
            default_flow_style=False,
        )


######################
## Model Criticism
######################


def get_model_fit(result, criterion="AIC", nll=None, n=None, M=None):
    """
    Given a fit result from scipy.optimize.minimize, compute AIC or BIC.

    Args:
    -----
    n: int, number of data points
    nll: float, negative log-likelihood value
    """
    assert M is not None, "Must provide the number of components M"
    assert n is not None, "Must provide the number of data points n"
    if nll is None:
        nll = result.fun  # this is the negative log-likelihood
    # Count how many free parameters are in the model
    # from the code, it’s 3 (w, mu, sigma) + 2M (b_j, r_j)
    n_params = 3 + 3 * M
    # Compute AIC or BIC
    # AIC = 2 * nll + 2 * n_params
    # BIC = 2 * nll + n_params * ln(n)
    if criterion.upper() == "AIC":
        ic_value = 2.0 * nll + 2.0 * n_params
    elif criterion.upper() == "BIC":
        ic_value = 2.0 * nll + n_params * np.log(n)
    else:
        raise ValueError(f"Unsupported criterion: {criterion}. Must be 'AIC' or 'BIC'.")
    fit_info = {
        "M": M,
        "nll": nll,
        "n_params": n_params,
        "IC": ic_value,
        "result": result,
    }
    return fit_info


def pick_best_model(max_M, pick_criterion, logs, holdoff_rate):
    assert pick_criterion in [
        "nll",
    ], f"Unknown criterion {pick_criterion}"
    suffix = "test" if holdoff_rate > 0 else "train"
    score_name = f"{pick_criterion.lower()}_{suffix}"
    best_M = None
    best_W = None
    for M in range(1, max_M + 1):
        # if logs[M] does not exist, skip
        if M not in logs:
            continue
        if best_M is None or logs[M][score_name] < logs[best_M][score_name]:
            best_M = M
            best_W = logs[M]["w_est"]
    return best_M, best_W




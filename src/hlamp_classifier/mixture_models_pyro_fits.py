# Given a fitted model, extract X, Y coordinates for plotting

import numpy as np
import pandas as pd
import yaml
import torch
import pyro.distributions as dist
import os
import argparse

from mixture_models_pyro_model import shifted_lognormal_pdf, gaussian_mixture_pdf
from mixture_models_pyro_utils import load_data
from mixture_models_pyro_globals import GLOBAL_GET_MAX_GAUSS_SIGMA


def extract_fit(best_model_path, conf_path, num_points=100, x_range=None):
    """
    Given a fitted model, extract

    Parameters
    ----------
    best_model_path : str
        Path to the best model. A yaml file with the best model parameters.
    conf_path : str
        Path to the yaml file with the configs, including the data path.
    num_points : int
        Number of points to use in the x_range.
    x_range : tuple
        The range of x values to use.

    Returns
    -------
    x_range : torch.Tensor
        The x values to plot.
    full_pdf : torch.Tensor
        The full pdf.
    weighted_ln : torch.Tensor
        The weighted lognormal pdf.
    weighted_normals : torch.Tensor
        The weighted normal pdf.
    """
    # Load the yaml configs for each
    with open(conf_path, "r") as f:
        configs = yaml.load(f, Loader=yaml.FullLoader)

    with open(best_model_path, "rb") as f:
        best_model_config = yaml.load(f, Loader=yaml.FullLoader)

    effective_width = configs["effective_width"]

    MAX_GAUSS_SIGMA = GLOBAL_GET_MAX_GAUSS_SIGMA(effective_width)
    
    data, _, _ = load_data(configs["dat_path"], configs["locus"])
    
    w_est = None
    if 'w_est' in best_model_config:
        w_est = best_model_config['w_est'] 
        muLN_est = best_model_config['muLN_est']
        sigmaLN_est = best_model_config['sigmaLN_est']
        locLn_est = best_model_config['locLn_est']
    
    alphas_est = best_model_config['alphas_est']
    mus_est = best_model_config['mus_est']
    sigmas_est = best_model_config['sigmas_est']

    # Define a grid of x values.
    if x_range is None:
        x_min = data.min()
        x_max = data.max()
        x_range = torch.linspace(x_min, x_max, num_points)
    else:
        x_range = torch.linspace(x_range[0], x_range[1], num_points)
    pdf_normals = gaussian_mixture_pdf(x_range, torch.tensor(alphas_est), torch.tensor(mus_est), torch.tensor(sigmas_est))

    if w_est is None:
        full_pdf = None
        weighted_ln = None
        weighted_normals = pdf_normals
    else:
        pdf_ln = shifted_lognormal_pdf(x_range, muLN_est, sigmaLN_est, locLn_est)
        weighted_ln = w_est * pdf_ln
        weighted_normals = (1 - w_est) * pdf_normals
        full_pdf = weighted_ln + weighted_normals

    return x_range, full_pdf, weighted_ln, weighted_normals


def export_model_fits(best_model_path, conf_path, num_points=100, tag=None, outDir=None):
    """
    Export the model fits, i.e., weighted_normals and the weighted_ln to a csv file.
    
    Parameters
    ----------
    best_model_path : str
        Path to the best model. A yaml file with the best model parameters.
    conf_path : str
        Path to the yaml file with the configs, including the data path.
    num_points : int
        Number of points to use in the x_range.
    tag : str
        Tag to add to the output file name.
    outDir : str
        Output directory. If None, use the directory of the best model path.
    
    Returns
    -------
    None
    """
    fname = 'evals.csv' if tag is None else f'evals_{tag}.csv'
    fpath = os.path.join(outDir, 'evals', fname)
    os.makedirs(os.path.dirname(fpath), exist_ok=True)

    x_range, full_pdf, weighted_ln, weighted_normals = extract_fit(
        best_model_path, conf_path, num_points
    )
    if outDir is None:
        outDir = os.path.dirname(best_model_path)

    # Save the weighted_normals and the weighted_ln as csv
    df_evals = pd.DataFrame(
        {"x": x_range, "weighted_ln": weighted_ln, "weighted_normals": weighted_normals}
    )
    df_evals.to_csv(fpath, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--model_path", type=str, help="Path to the model output. Ignored if best_model_path or conf_path are provided", required=False)
    parser.add_argument("--best_model_path", type=str, help="Path to the best model", required=False)
    parser.add_argument("--conf_path", type=str, help="Path to the configs file", required=False)
    parser.add_argument(
        "--num_points",
        type=int,
        default=100,
        help="Number of points to use in the x_range",
    )
    parser.add_argument("--outDir", type=str, default=None, help="Output directory")
    args = parser.parse_args()

    # Assert that either model_path or best_model_path and conf_path are provided
    if args.model_path is not None:
        assert args.best_model_path is None and args.conf_path is None, "If model_path is provided, best_model_path and conf_path should not be provided"
        args.best_model_path = os.path.join(args.model_path, "best_model.yaml")
        args.conf_path = os.path.join(args.model_path, "configs.yaml")
    elif args.best_model_path is None or args.conf_path is None:
        raise ValueError("Either model_path or best_model_path and conf_path should be provided")
    print(f"Using best_model_path: {args.best_model_path}")
    print(f"Using conf_path: {args.conf_path}")
    export_model_fits(
        best_model_path=args.best_model_path,
        conf_path=args.conf_path,
        num_points=args.num_points,
        outDir=args.outDir,
    )

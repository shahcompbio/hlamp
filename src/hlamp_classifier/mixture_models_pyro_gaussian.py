import pickle
import os
import time
import argparse
import yaml
from tqdm import tqdm
import torch
import pyro
import pyro.distributions as dist
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyro.optim as optim
from pyro.infer import SVI, Trace_ELBO, Predictive
from mixture_models_pyro_utils import load_data
from mixture_models_pyro_fits import export_model_fits
from mixture_utils import write_results_to_yaml, pick_best_model
from mixture_models_pyro_model import (
    make_init_to_custom_value,
    model_gaussian,
    gaussian_mixture_pdf,
    model,
)  # now model is the new Gaussian mixture model
from mixture_models_pyro_globals import (
    GLOBAL_GET_MIN_SEP,
    GLOBAL_GET_MAX_GAUSS_SIGMA,
)
from mixture_models_pyro_globals import get_effective_width

# ------------------------------------------------------------------------------
# Set device to GPU if available
# ------------------------------------------------------------------------------
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if device.type == "cuda":
    torch.set_default_tensor_type(torch.cuda.FloatTensor)
    print("Using GPU")
else:
    print("Not using GPU, using CPU")


# ------------------------------------------------------------------------------
# Helper: Get parameter value either from the param store or from a provided guide sample.
# ------------------------------------------------------------------------------
def get_param(name, local_store=None, guide_sample=None):
    if guide_sample is not None:
        return guide_sample[name]
    if local_store is None:
        ps = pyro.get_param_store()
        key = f"AutoDelta.{name}"
        if key in ps:
            return ps[key]
        elif name in ps:
            return ps[name]
        else:
            raise KeyError(
                f"Parameter {name} not found in param store. Available keys: {list(ps.keys())}"
            )
    else:
        return local_store[name]


# ------------------------------------------------------------------------------
# Compute negative log-likelihood using the inferred Gaussian mixture parameters.
# ------------------------------------------------------------------------------
def compute_nll(data, MAX_GAUSS_SIGMA=None, effective_width=None, fit_result=None):
    """
    Constructs the Gaussian mixture using the inferred parameters and computes the NLL.
    """
    # Unpack parameters: alpha, ordered means (mus) and standard deviations (sigmas)
    alphas_est, mus_est, sigmas_est, delta_est = fit_result.unpack_params(
        keep_tensor=True,
    )
    # Build the Gaussian components list.
    M = len(alphas_est)
    gauss_components = [dist.Normal(mus_est[i], sigmas_est[i]) for i in range(M)]
    mix = HeteroMixture(alphas_est, gauss_components)
    data_tensor = torch.tensor(data, dtype=torch.float)
    nll = -mix.log_prob(data_tensor).sum()
    return nll.detach().item()


# ------------------------------------------------------------------------------
# A function to plot the raw data and the fitted Gaussian mixture density.
# ------------------------------------------------------------------------------
def internal_plot(data, x_range, mixture_pdf, extra_pdfs=None, fpath=None, titleStr=None):
    plt.clf()
    plt.figure(figsize=(10, 6))
    plt.hist(data.numpy(), bins=30, density=True, alpha=0.5, label="Data")
    if x_range is not None:
        if mixture_pdf is not None:
            plt.plot(
                x_range.numpy(),
                mixture_pdf.detach().numpy(),
                color="blue",
                label="Gaussian Mixture",
            )
        # If extra_pdfs is provided, plot them as well.
        if extra_pdfs is not None:
            extra_colors = ["red", "green", "orange", "purple"]
            for key, pdf in extra_pdfs.items():
                plt.plot(
                    x_range.numpy(),
                    pdf.detach().numpy(),
                    color=extra_colors.pop(0),
                    linestyle="--",
                    label=key,
                )
    # Make the labels of x axes the original, not the log values
    if data.max() < 10:
        xticks = plt.xticks()[0]
        xtick_labels = [f"{np.exp(x):.1f}" for x in xticks]
        plt.xticks(xticks, labels=xtick_labels)
    plt.xlabel("Copy Number")
    plt.ylabel("Density")
    plt.legend()
    min_title = f"Fitted Mixture Model\n{len(data)} data points"
    full_title = f"{titleStr}\n{min_title}" if titleStr else min_title
    plt.title(full_title)
    if fpath is None:
        plt.show()
    else:
        plt.savefig(fpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_fit(
    data,
    x_range=None,
    num_points=500,
    MAX_GAUSS_SIGMA=None,
    fpath=None,
    effective_width=None,
    fit_result=None,
    full_data=None,
    titleStr=None,
):
    """
    Plots the histogram of the data and the fitted Gaussian mixture density.
    """
    # Retrieve fitted parameters.
    alphas_est, mus_est, sigmas_est, delta_est = fit_result.unpack_params(
        keep_tensor=True,
    )
    # Print model fits
    print("##############################################")
    print("Model Fits:")
    print("##############################################")
    alpha_print = (
        alphas_est.detach().cpu().numpy()
        if hasattr(alphas_est, "detach")
        else alphas_est
    )
    mu_print = mus_est.detach().cpu().numpy() if hasattr(mus_est, "detach") else mus_est
    sigma_print = (
        sigmas_est.detach().cpu().numpy()
        if hasattr(sigmas_est, "detach")
        else sigmas_est
    )
    delta_print = (
        delta_est.detach().cpu().numpy() if hasattr(delta_est, "detach") else delta_est
    )
    print(f"alphas_est: {np.round(alpha_print, 2)}")
    print(f"mus_est: {np.round(mu_print, 2)}")
    print(f"sigmas_est: {np.round(sigma_print, 2)}")
    print(f"delta_est: {np.round(delta_print, 2)}")
    print(f"max cn: {data.max():.2f}")
    print(f"sigma0: {fit_result.predictive['sigma0']:.2f}")
    print("##############################################")
    # Define a grid of x values.
    if x_range is None:
        x_min = data.min()
        x_max = data.max()
        x_range = torch.linspace(x_min, x_max, num_points)
    else:
        x_range = torch.linspace(x_range[0], x_range[1], num_points)
    # Compute the overall Gaussian mixture density.
    # mixture_pdf = gaussian_mixture_pdf(x_range, alphas_est, mus_est, sigmas_est)
    M = len(alphas_est)
    gauss_components = [dist.Normal(mus_est[i], sigmas_est[i]) for i in range(M)]
    mix = HeteroMixture(alphas_est, gauss_components)
    mixture_pdf = mix.log_prob(x_range).detach().exp()
    # mix = dist.MixtureSameFamily(
    #     dist.Categorical(alphas_est),
    #     dist.Normal(mus_est, sigmas_est),
    # )
    # Now compute the individual Gaussian components.
    extra_pdfs = {}
    for i in range(M):
        # gauss_pdf = dist.Normal(mus_est[i], sigmas_est[i]).log_prob(x_range).detach().exp()
        # mix_ = HeteroMixture([alphas_est[i]], [gauss_components[i]])
        # gauss_pdf = mix_.log_prob(x_range).detach().exp()
        gauss_pdf = alphas_est[i] * gauss_components[i].log_prob(x_range).detach().exp()
        extra_pdfs[f"gauss_comp_{i}"] = gauss_pdf

    # Now fit a Gaussian to the data using optimization and sklearn
    from scipy.stats import norm

    mu, sd = norm.fit(data.numpy())
    # Print the fitted parameters
    print(f"Fitted Single Gaussian: mu = {mu:.2f}, sd = {sd:.2f}")
    # Create a normal distribution with the fitted parameters
    fitted_dist = dist.Normal(mu, sd)
    # Compute the PDF of the fitted distribution
    fitted_pdf = fitted_dist.log_prob(x_range).detach().exp()
    # Append the fitted PDF to the extra_pdfs list
    extra_pdfs["Single_Gaussian"] = fitted_pdf
    internal_plot(
        data=data,
        x_range=x_range,
        mixture_pdf=mixture_pdf,
        extra_pdfs=extra_pdfs,
        fpath=fpath,
        titleStr=titleStr,
    )
    # Just plot the non log version too
    fpath_orig = fpath.replace(".png", "_orig.png")
    internal_plot(
        data=np.exp(data),
        x_range=None,
        mixture_pdf=None,
        extra_pdfs=extra_pdfs,
        fpath=fpath_orig,
        titleStr=titleStr,
    )
    if full_data is not None:
        fpath_orig_full = fpath.replace(".png", "_orig_full.png")
        internal_plot(
            data=np.exp(full_data),
            x_range=None,
            mixture_pdf=None,
            extra_pdfs=extra_pdfs,
            fpath=fpath_orig_full,
            titleStr=titleStr,
        )

    def fit_gmm(data=None, covariance_type="diag"):
        """Fit a mixture of gaussians with two components and print the means and sigmas and the components"""
        from sklearn.mixture import GaussianMixture

        # gmm = GaussianMixture(n_components=2, covariance_type=covariance_type, max_iter=10000, n_init=1, random_state=111)
        gmm = GaussianMixture(n_components=2, covariance_type=covariance_type)
        gmm.fit(data.reshape(-1, 1))
        # Get the components
        means = np.squeeze(gmm.means_)
        sds = np.squeeze(np.sqrt(gmm.covariances_))
        weights = np.squeeze(gmm.weights_)
        print("GMM fit is ")
        print(
            f"Means = {np.round(means, 2)} -- sds = {np.round(sds, 2)} -- weights = {np.round(weights, 2)}"
        )
        return means, sds, weights

    # Fit a GMM to the data
    means, sds, weights = fit_gmm(data, covariance_type="diag")
    # Plot them
    extra_pdfs = {}
    try:
        gauss_components = [dist.Normal(means[i], sds[i]) for i in range(M)]
        for i in range(M):
            gauss_pdf = weights[i] * gauss_components[i].log_prob(x_range).exp()
            extra_pdfs[f"GMM_comp_{i}"] = gauss_pdf
    except Exception as e:
        print(f"Failed to compute GMM components: {e}")
    fpath_gmm = fpath.replace(".png", "_gmm.png")
    internal_plot(
        data=data,
        x_range=x_range,
        mixture_pdf=None,
        extra_pdfs=extra_pdfs,
        fpath=fpath_gmm,
        titleStr=titleStr,
    )


# ------------------------------------------------------------------------------
# A custom heterogeneous mixture distribution (unchanged).
# ------------------------------------------------------------------------------
class HeteroMixture(dist.TorchDistribution):
    arg_constraints = {}
    support = dist.constraints.real
    has_rsample = False

    def __init__(self, mixing_probs, components, validate_args=None):
        self.mixing_probs = mixing_probs
        self.components = components
        batch_shape = components[0].batch_shape
        event_shape = components[0].event_shape
        super(HeteroMixture, self).__init__(
            batch_shape=batch_shape,
            event_shape=event_shape,
            validate_args=validate_args,
        )

    def log_prob(self, value):
        log_probs = []
        for i, comp in enumerate(self.components):
            log_w = torch.log(self.mixing_probs[i] + 1e-8)
            log_probs.append(log_w + comp.log_prob(value))
        log_probs = torch.stack(log_probs, dim=0)
        return torch.logsumexp(log_probs, dim=0)


# ------------------------------------------------------------------------------
# The FitResult class now unpacks only the Gaussian mixture parameters.
# ------------------------------------------------------------------------------
class FitResult:
    def __init__(
        self,
        success,
        message,
        x,
        fun,
        MAX_GAUSS_SIGMA=None,
        effective_width=None,
        predictive=None,
        guide_sample=None,
    ):
        self.success = success
        self.message = message
        self.x = x
        self.fun = fun
        self.MAX_GAUSS_SIGMA = MAX_GAUSS_SIGMA
        self.effective_width = effective_width
        self.predictive = predictive
        self.guide_sample = guide_sample
        self.set_results()

    def set_results(self):
        alphas_est, mus_est, sigmas_est, delta_est = self.unpack_params()
        self.alpha = alphas_est
        self.mus_est = mus_est
        self.sigmas_est = sigmas_est
        self.delta_est = delta_est

    def unpack_params(
        self,
        keep_tensor=False,
    ):
        """ """
        if self.predictive is not None:
            alphas_est = self.predictive["alphas"]
            mus_est = self.predictive["mus_est"]
            sigmas_est = self.predictive["sigma_est"]
            delta_est = self.predictive["delta"]
        elif self.guide_sample is not None:
            alphas_est = self.guide_sample["alphas_est"]
            mus_est = self.guide_sample["mus_est"]
            sigmas_est = self.guide_sample["sigma_est"]
            delta_est = self.guide_sample["delta"]
        else:
            raise ValueError(
                "Either predictive or guide_sample must be provided to unpack_params."
            )
        if not keep_tensor:
            alphas_est = alphas_est.detach().cpu().numpy()
            mus_est = mus_est.detach().cpu().numpy()
            sigmas_est = sigmas_est.detach().cpu().numpy()
            delta_est = delta_est.detach().cpu().numpy()

        return alphas_est, mus_est, sigmas_est, delta_est


# ------------------------------------------------------------------------------
# The fit_model function updated to work with the new Gaussian mixture model.
# ------------------------------------------------------------------------------
def fit_model(
    M,
    train_counts,
    test_counts,
    outDir=None,
    figuresDir=None,
    titleStr=None,
    bfb_family="gaussian",  # updated family name
    verbose=False,
    effective_width=None,
    n_steps=1000,
    learning_rate=0.005,
    tag="",
    seed=42,
    use_init=True,
    r_counts=None,
):
    MAX_GAUSS_SIGMA = GLOBAL_GET_MAX_GAUSS_SIGMA(effective_width)
    print(f"Number of cells in the train set: {len(train_counts)}")
    print(
        f"Number of cells in the test set: {len(test_counts) if test_counts is not None else 0}"
    )
    print(f"Effective width: {effective_width:.2f}")
    print(f"MAX_GAUSS_SIGMA: {MAX_GAUSS_SIGMA:.2f}")
    print(f"Data mean: {train_counts.mean():.2f}")
    print(f"Data std: {train_counts.std():.2f}")
    print(f"Data kurtosis: {pd.Series(np.exp(train_counts)).kurtosis():.2f}")
    pyro.clear_param_store()
    pyro.set_rng_seed(seed)

    def make_init_func(train_data, M_value):
        def init_to_custom_value(site):
            name = site["name"]
            # TODO: fix to alphas
            if name == "alpha":
                # return torch.ones(M_value) / M_value
                return torch.ones(M_value) * 0.2
            elif name == "mu_components":
                return torch.ones(M_value) * train_data.mean()
            # elif name == "l_components":
            #     return torch.zeros(M_value)
            else:
                # For any other latent variables, fallback to the default initialization
                return site["fn"]()

        return init_to_custom_value

    def make_init_xx(train_data, M_value, effective_width=None):
        """
        Custom initialization function for the Gaussian mixture model.

        Parameters:
        train_data: the training data, used to compute the data mean.
        M_value: number of Gaussian components.
        effective_width: (optional) if provided, you could use this to adjust the baseline separation.

        Returns:
        A function that maps a site (dictionary) to an initial value.
        """
        # You can replace this constant with a call to GLOBAL_GET_MIN_SEP(effective_width) if desired.
        baseline = 1.5

        def init_to_custom_value(site):
            name = site["name"]
            if name == "alphas":
                # Initialize mixing weights uniformly.
                return torch.ones(M_value) / M_value
                # return torch.tensor(np.array([0.01, .09]))
            # elif name == "mu0":
            #     # Initialize the first mean near the data mean (with small noise).
            #     return torch.tensor(train_data.mean()) + 0.1 * torch.randn(1)
            # elif name == "delta":
            #     # Initialize M-1 positive differences. The prior expectation of an Exponential(1) is 1.
            #     # These will be shifted by 'baseline' in the model.
            #     # return 1.0 * torch.ones(M_value - 1) + 0.1 * torch.rand(M_value - 1)
            #     return 0.1 * torch.rand(M_value - 1)
            # elif name == "l_components":
            #     # Initialize near zero so that sigma_components starts at ~0.5*MAX_GAUSS_SIGMA.
            #     return torch.zeros(M_value)
            else:
                # Fallback to the default initialization.
                return site["fn"]()

        return init_to_custom_value

    # init_loc_fn = make_init_func(train_counts, M)
    # guide = pyro.infer.autoguide.AutoDelta(model_gaussian, init_loc_fn=init_loc_fn)
    # guide = pyro.infer.autoguide.AutoDelta(model_gaussian, init_loc_fn=make_init_xx(train_counts, M))
    guide = pyro.infer.autoguide.AutoDelta(model_gaussian)
    optimizer = optim.AdamW({"lr": learning_rate})
    svi = SVI(model_gaussian, guide, optimizer, loss=Trace_ELBO())
    # guide = pyro.infer.autoguide.AutoDelta(model)
    # optimizer = optim.AdamW({"lr": learning_rate})
    # svi = SVI(model, guide, optimizer, loss=Trace_ELBO())

    best_loss = float("inf")
    best_state = None
    worst_loss = float("-inf")
    pbar = tqdm(range(n_steps))
    for step in pbar:
        loss = svi.step(
            train_counts,
            MAX_GAUSS_SIGMA=MAX_GAUSS_SIGMA,
            M=M,
            effective_width=effective_width,
        )
        if step == 0:
            init_loss = loss
        if loss > worst_loss:
            worst_loss = loss
        if loss < best_loss:
            best_loss = loss
            best_state = pyro.get_param_store().get_state()
        pbar.set_description(
            f"Step {step} : loss = {loss:.2f} -- best_loss = {best_loss:.2f} -- worst_loss = {worst_loss:.2f} -- init_loss = {init_loss:.2f}"
        )

    if loss > best_loss and best_state is not None:
        pyro.get_param_store().set_state(best_state)
        print(f"Restored the best model with loss = {best_loss:.2f}")

    # Print all parameters in the param store.
    predictive_svi = Predictive(model_gaussian, guide=guide, num_samples=1)(
        train_counts,
        MAX_GAUSS_SIGMA=MAX_GAUSS_SIGMA,
        M=M,
        effective_width=effective_width,
    )

    # Sqeeuze everything in predictive
    new_pred = {}
    for k, v in predictive_svi.items():
        if isinstance(v, torch.Tensor):
            new_pred[k] = torch.squeeze(v)
        else:
            new_pred[k] = v
    fit_result = FitResult(
        success=True,
        message="Success",
        x=train_counts,
        fun=None,
        MAX_GAUSS_SIGMA=MAX_GAUSS_SIGMA,
        effective_width=effective_width,
        predictive=new_pred,
        guide_sample=None,
    )

    nll = compute_nll(
        train_counts,
        MAX_GAUSS_SIGMA=MAX_GAUSS_SIGMA,
        effective_width=effective_width,
        fit_result=fit_result,
    )

    fit_result.fun = nll

    fPath = os.path.join(figuresDir, f"fit_result_M{M}{tag}.png")
    print(f"alphas_est: {fit_result.alpha}")
    nll_train = fit_result.fun

    #titleStr = f"{titleStr}\nMixture weights: {np.round(fit_result.alpha, 2)}"
    # Concatenate train and test count as a torch.tensor
    full_data = None
    if test_counts is not None:
        full_data = r_counts
        full_data = torch.tensor(full_data, dtype=torch.float)
    plot_fit(
        torch.tensor(train_counts, dtype=torch.float),
        fpath=fPath,
        MAX_GAUSS_SIGMA=MAX_GAUSS_SIGMA,
        effective_width=effective_width,
        fit_result=fit_result,
        full_data=full_data,
        titleStr=titleStr,
    )

    tmp_log = {"M": M}
    tmp_log["alphas_est"] = fit_result.alpha
    tmp_log["mus_est"] = fit_result.mus_est
    tmp_log["sigmas_est"] = fit_result.sigmas_est
    tmp_log["nll_train"] = nll_train
    tmp_log["success"] = fit_result.success
    tmp_log["message"] = fit_result.message
    tmp_log["loss"] = loss

    if test_counts is not None:
        nll_test = compute_nll(
            test_counts,
            MAX_GAUSS_SIGMA=MAX_GAUSS_SIGMA,
            effective_width=effective_width,
            fit_result=fit_result,
        )
    else:
        nll_test = None
    tmp_log["nll_test"] = nll_test

    return tmp_log, fit_result




# ------------------------------------------------------------------------------
# Finding the smallest x that covers 95% of the data (kinda the range based on alpha)
# ------------------------------------------------------------------------------

def component_quantiles(fit_result, p=0.95):
    mu=torch.as_tensor(fit_result.mus_est,dtype=torch.float)
    sg=torch.as_tensor(fit_result.sigmas_est,dtype=torch.float)
    z=Normal(0.0,1.0).icdf(torch.tensor(p,dtype=torch.float))
    return (mu+sg*z)

def mixture_quantile(fit_result, p=0.95, data=None, pad=2.0, iters=60):
    w=torch.as_tensor(fit_result.alpha,dtype=torch.float)
    mu=torch.as_tensor(fit_result.mus_est,dtype=torch.float)
    sg=torch.as_tensor(fit_result.sigmas_est,dtype=torch.float)
    Phi=Normal(0.0,1.0).cdf
    if data is not None:
        lo=float(torch.min(data))-pad*float(torch.std(data))
        hi=float(torch.max(data))+pad*float(torch.std(data))
    else:
        lo=float(torch.min(mu-6*sg))
        hi=float(torch.max(mu+6*sg))
    lo=torch.tensor(lo); hi=torch.tensor(hi)
    for _ in range(iters):
        mid=(lo+hi)/2
        cdf=torch.sum(w*(Phi((mid-mu)/sg)))
        lo,hi=torch.where(cdf<p,mid,lo),torch.where(cdf>=p,mid,hi)
    return float((lo+hi)/2)


# ------------------------------------------------------------------------------
# Finding the mode
# ------------------------------------------------------------------------------
import torch
import math
from torch.distributions import Normal, Categorical

import torch
from torch.distributions import Normal

def mixture_logpdf(x, weights, mus, sigmas):
    # x: scalar tensor or 1D; returns scalar log-density at x
    comp = Normal(mus, sigmas)
    log_terms = torch.log(weights) + comp.log_prob(x.unsqueeze(-1))
    return torch.logsumexp(log_terms, dim=-1)

def find_mode_1d(fit_result, data, pad=1.0, n_grid=2048, n_refine=80, lr=0.1, device=None):
    # gather tensors on a common device
    device = device or (torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu"))
    w  = torch.as_tensor(fit_result.alpha,     dtype=torch.float, device=device)
    mu = torch.as_tensor(fit_result.mus_est,   dtype=torch.float, device=device)
    sg = torch.as_tensor(fit_result.sigmas_est,dtype=torch.float, device=device)
    data = data.to(device)
    # coarse grid init on logpdf
    lo = data.min().item() - pad * data.std().item()
    hi = data.max().item() + pad * data.std().item()
    grid = torch.linspace(lo, hi, n_grid, device=device)
    logvals = mixture_logpdf(grid, w, mu, sg)
    x0 = grid[torch.argmax(logvals)].detach()

    # refine with gradient ascent on log-density
    x = torch.tensor(x0.item(), device=device, requires_grad=True)
    opt = torch.optim.Adam([x], lr=lr)
    for _ in range(n_refine):
        opt.zero_grad()
        loss = -mixture_logpdf(x, w, mu, sg)  # maximize logpdf
        loss.backward()
        opt.step()
        # keep x in a sane region
        with torch.no_grad():
            x.clamp_(lo, hi)
    x_mode = x.detach()
    f_xstar = torch.exp(mixture_logpdf(x_mode, w, mu, sg)).item()
    return x_mode, f_xstar

@torch.no_grad()
def mass_in_window_around_mode(fit_result, x_mode, r):
    w  = torch.as_tensor(fit_result.alpha,     dtype=torch.float)
    mu = torch.as_tensor(fit_result.mus_est,   dtype=torch.float)
    sg = torch.as_tensor(fit_result.sigmas_est,dtype=torch.float)
    Phi = Normal(0.0, 1.0).cdf
    z_hi = (x_mode + r - mu) / sg
    z_lo = (x_mode - r - mu) / sg
    return torch.sum(w * (Phi(z_hi) - Phi(z_lo))).item()


@torch.no_grad()
def responsibilities_at_x(fit_result, x):
    w = torch.as_tensor(fit_result.alpha, dtype=torch.float)
    mu = torch.as_tensor(fit_result.mus_est, dtype=torch.float)
    sg = torch.as_tensor(fit_result.sigmas_est, dtype=torch.float)
    comp = Normal(mu, sg)
    log_w_post = torch.log(w) + comp.log_prob(torch.as_tensor(x))
    log_w_post = log_w_post - torch.logsumexp(log_w_post, dim=-1)
    return torch.softmax(log_w_post, dim=-1).numpy()

# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------
def main(
    dat_path=None,
    locus=None,
    subset_counts=0,
    seed=42,
    max_M=5,
    holdoff_rate=0.2,
    outDir=None,
    datName=None,
    pick_criterion="nll",
    no_loop=False,
    bfb_family="gaussian",
    verbose=False,
    n_steps=1000,
    learning_rate=0.005,
):
    LOG_PATH = os.path.join(outDir, "logs.yaml")
    TRAINING_LOG_PATH = os.path.join(outDir, "training_logs.yaml")
    BEST_MODEL_PATH = os.path.join(outDir, "best_model.yaml")
    CONF_PATH = os.path.join(outDir, "configs.yaml")
    BEST_MODEL_PATH_FULL = os.path.join(outDir, "best_models", "best_model_full.yaml")
    os.makedirs(os.path.dirname(BEST_MODEL_PATH_FULL), exist_ok=True)

    start_time = time.time()
    outDir = os.path.abspath(outDir)
    figuresDir = os.path.join(outDir, "figures")
    os.makedirs(figuresDir, exist_ok=True)

    if datName is None:
        datName = os.path.basename(dat_path)
        datName = os.path.splitext(datName)[0]

    # Remove anything after . from the datName
    datNameTitle = datName.split(".")[0]
    titleStr = f"{datNameTitle} {locus}"
    print(f"Running for {titleStr}")

    configs = {
        "dat_path": dat_path,
        "locus": locus,
        "subset_counts": subset_counts,
        "seed": seed,
        "max_M": max_M,
        "holdoff_rate": holdoff_rate,
        "outDir": outDir,
        "datName": datName,
        "pick_criterion": pick_criterion,
        "no_loop": no_loop,
        "bfb_family": bfb_family,
        "n_steps": n_steps,
    }
    with open(CONF_PATH, "w") as f:
        yaml.dump(configs, f)
    
    logs = {}
    r_counts, tail_vals, r_counts_FULL = load_data(dat_path, locus=locus, do_log=True)

    n_cells_orig = len(r_counts)
    # If there are fewer than 10 cells don't run anything just exit
    if n_cells_orig < 10:
        print(
            f"Not enough cells ({n_cells_orig}) to run the model. Exiting. Please check your data."
        )
        return
    np.random.seed(seed)
    torch.manual_seed(seed)
    print(f"Loaded {len(r_counts)} cells")
    configs["n_cells"] = n_cells_orig

    if n_cells_orig < 150 and holdoff_rate > 0.05:
        print(
            f"Holdoff rate {holdoff_rate} is too high for {n_cells_orig} cells. Setting to 0.05"
        )
        holdoff_rate = 0.05
    if subset_counts > 0:
        print(f"Subsetting to {subset_counts} cells")
        n_cells = len(r_counts)
        sel_idx = np.random.choice(
            n_cells, size=min(subset_counts, n_cells), replace=False
        )
        r_counts = r_counts[sel_idx]

    if holdoff_rate > 0:
        assert holdoff_rate < 1, "Holdoff rate should be between 0 and 1"
        n = len(r_counts)
        n_train = int(n * (1 - holdoff_rate))
        np.random.shuffle(r_counts)
        train_counts = r_counts[:n_train]
        test_counts = r_counts[n_train:]
    else:
        train_counts = r_counts
        test_counts = None

    if tail_vals is not None:
        train_counts = np.concatenate([train_counts, tail_vals])

    effective_width = get_effective_width(train_counts)
    configs["effective_width"] = float(effective_width)


    # Plot r_counts_full
    if r_counts_FULL is not None:
        fpath = os.path.join(figuresDir, f"r_counts_FULL_{datNameTitle}_{locus}.png")
        internal_plot(
            data=torch.tensor(r_counts_FULL, dtype=torch.float),
            x_range=None,
            mixture_pdf=None,
            extra_pdfs=None,
            fpath=fpath,
            titleStr=titleStr,
        )

    if no_loop:
        M = max_M
        print(f"Fitting model with M = {M}...")
        losses = {}
        tmp_logs = {}
        the_fits = {}
        for the_seed in [10, 20, 30, 40, 50]:
        #for the_seed in [50]:
            # for the_seed in [seed]:
            print(f"Running for seed {the_seed}")
            # set the seed for all
            np.random.seed(the_seed)
            torch.manual_seed(the_seed)
            tmp_log, fit_result = fit_model(
                M=M,
                outDir=outDir,
                train_counts=train_counts,
                test_counts=test_counts,
                figuresDir=figuresDir,
                titleStr=titleStr + f" M = {M}",
                bfb_family=bfb_family,
                verbose=verbose,
                effective_width=effective_width,
                n_steps=n_steps,
                learning_rate=learning_rate,
                seed=the_seed,
                use_init=False,
                r_counts=r_counts,
            )
            ##############################################################################
            # Extract the mode
            ##############################################################################
            data_t = torch.as_tensor(train_counts, dtype=torch.float)
            x_star, f_xstar = find_mode_1d(fit_result, data_t)  # now works
            p_window = mass_in_window_around_mode(fit_result, x_star.item(), r=0.5)
            resp = responsibilities_at_x(fit_result, x_star.item())
            print({"x_mode": float(x_star), "density_at_mode": f_xstar, "mass_in_window": p_window, "resp": resp})
            ##############################################################################
            # Extract the model range
            ##############################################################################
            q_comp = component_quantiles(fit_result, p=0.95)  # tensor of size M
            # overall mixture 95% upper coverage point
            q_mix = mixture_quantile(fit_result, p=0.95, data=torch.as_tensor(train_counts,dtype=torch.float))
            print(f"Component 95% quantiles: {q_comp}")
            print(f"Overall mixture 95% quantile: {q_mix:.2f}")
            # Add them to tmp_log
            tmp_log["x_mode"] = float(x_star)
            tmp_log["density_at_mode"] = float(f_xstar)
            tmp_log["mass_in_window_0.5"] = float(p_window)
            tmp_log["resp"] = resp
            tmp_log["comp_95_quantiles"] = q_comp.numpy().tolist()
            tmp_log["mix_95_quantile"] = float(q_mix)
            tmp_logs[the_seed] = tmp_log
            losses[the_seed] = tmp_log["loss"]
            #losses[the_seed] = tmp_log["nll_test"]
            the_fits[the_seed] = fit_result
        
        # data_t = torch.as_tensor(train_counts, dtype=torch.float)
        # x_star, f_xstar = find_mode_1d(fit_result, data_t)    # x_star is the mode, f_xstar is the density there
        # r = 0.5                                               # choose neighborhood radius in data units
        # p_window = mass_in_window_around_mode(fit_result, x_star.item(), r)
        #resp = responsibilities_at_x(fit_result, x_star.item())
        #print({"x_mode": float(x_star), "density_at_mode": f_xstar, "mass_in_window": p_window, "resp": resp})        

        best_seed = min(losses, key=losses.get)
        logs[M] = tmp_logs[best_seed]
        best_fit = the_fits[best_seed]
        write_results_to_yaml(tmp_logs, TRAINING_LOG_PATH, verbose=verbose)
        # Print the chosen alpha
        print(f"Best alpha_est = {logs[M]['alphas_est']}")
        # plot the latest
        #titleStr = f"{titleStr}\nMixture weights: {np.round(fit_result.alpha, 2)}"
        plot_fit(
            torch.tensor(train_counts, dtype=torch.float),
            fpath=os.path.join(outDir, f"fit_result_M{M}.png"),
            MAX_GAUSS_SIGMA=None,
            effective_width=effective_width,
            fit_result=best_fit,
            titleStr=titleStr,
        )
        try:
            write_results_to_yaml(logs, LOG_PATH)
            write_results_to_yaml(logs[M], os.path.join(outDir, "best_model.yaml"))
        except Exception as e:
            print(f"Failed to save logs: {e}")
    else:
        for M in range(1, max_M + 1):
            print(f"Fitting model with M = {M}...")
            try:
                losses = {}
                tmp_logs = {}
                for the_seed in [seed]:
                    np.random.seed(the_seed)
                    torch.manual_seed(the_seed)
                    tmp = fit_model(
                        M=M,
                        outDir=outDir,
                        train_counts=train_counts,
                        test_counts=test_counts,
                        figuresDir=figuresDir,
                        titleStr=titleStr + f" M = {M}",
                        bfb_family=bfb_family,
                        verbose=verbose,
                        effective_width=effective_width,
                        n_steps=n_steps,
                        learning_rate=learning_rate,
                        seed=the_seed,
                        use_init=True,
                        r_counts=r_counts,
                    )
                    tmp_logs[the_seed] = tmp
                    losses[the_seed] = tmp["loss"]
                best_seed = min(losses, key=losses.get)
                logs[M] = tmp_logs[best_seed]
            except Exception as e:
                print(f"Failed to fit model with M = {M}: {e}")
                continue
            try:
                write_results_to_yaml(logs, LOG_PATH)
            except Exception as e:
                print(f"Failed to save logs: {e}")
                continue

        best_M, best_W = pick_best_model(max_M, pick_criterion, logs, holdoff_rate)
        if best_M is not None and best_W is not None:
            print(f"Best model is M = {best_M} with weight = {best_W:.3f}")
            write_results_to_yaml(logs[best_M], os.path.join(outDir, "best_model.yaml"))
            print(logs[best_M])

    write_results_to_yaml(logs, LOG_PATH, verbose=verbose)
    print(f"Results are saved in {outDir}")
    end_time = time.time()
    print(f"Total time taken: {end_time - start_time:.2f} seconds")
    configs["total_time"] = end_time - start_time
    with open(CONF_PATH, "w") as f:
        yaml.dump(configs, f)

    export_model_fits(BEST_MODEL_PATH, CONF_PATH, outDir=outDir)

    # (Optional) Rerun on full data if needed...
    if M > 0 and holdoff_rate > 0 and not no_loop:
        best_model = logs[best_M]
        best_M = best_model["M"]
        print(f"Rerunning model with M = {best_M} on the full dataset...")
        losses = {}
        for new_seed in [100]:
            print(f"Running model with seed = {new_seed}...")
            np.random.seed(new_seed)
            torch.manual_seed(new_seed)
            np.random.shuffle(r_counts)
            tmp = fit_model(
                M=best_M,
                outDir=outDir,
                train_counts=r_counts,
                test_counts=None,
                figuresDir=figuresDir,
                titleStr=titleStr + f" M = {M} (full)",
                bfb_family=bfb_family,
                verbose=verbose,
                effective_width=effective_width,
                n_steps=n_steps,
                learning_rate=learning_rate,
                tag=f"_best_model_full_{new_seed}",
                seed=new_seed,
                use_init=True,
                r_counts=r_counts,
            )
            losses[new_seed] = tmp["loss"]
            BEST_MODEL_PATH_FULL_tmp = BEST_MODEL_PATH_FULL.replace(
                ".yaml", f"_{new_seed}.yaml"
            )
            write_results_to_yaml(tmp, BEST_MODEL_PATH_FULL_tmp)
            export_model_fits(
                BEST_MODEL_PATH_FULL_tmp,
                CONF_PATH,
                tag=f"full_{new_seed}",
                outDir=outDir,
            )
        best_seed = min(losses, key=losses.get)
        BEST_MODEL_PATH_FULL_tmp = BEST_MODEL_PATH_FULL.replace(
            ".yaml", f"_{best_seed}.yaml"
        )
        new_path = os.path.join(outDir, "best_model_full.yaml")
        os.system(f"cp {BEST_MODEL_PATH_FULL_tmp} {new_path}")
        os.system(
            f'cp {os.path.join(figuresDir, f"fit_result_M{best_M}_best_model_full_{best_seed}.png")} {os.path.join(outDir, "fit_result_best_model_full.png")}'
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the Gaussian mixture model on the data"
    )
    parser.add_argument("--dat_path", type=str, help="Path to the data file")
    parser.add_argument("--locus", type=str, help="Locus to analyze")
    parser.add_argument(
        "--subset_counts", type=int, default=0, help="Number of cells to subset"
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--max_M", type=int, default=6, help="Max number of components")
    parser.add_argument(
        "--holdoff_rate", type=float, default=0.2, help="Holdoff rate for test set"
    )
    parser.add_argument("--outDir", type=str, help="Output directory")
    parser.add_argument("--dataset_name", type=str, help="Name of the dataset")
    parser.add_argument(
        "--no_loop",
        action="store_true",
        help="Do not loop over M. Just run for the max M",
    )
    parser.add_argument(
        "--pick_criterion",
        type=str,
        default="nll",
        help="Criterion to pick the best model",
    )
    parser.add_argument(
        "--BFB_family",
        type=str,
        default="gaussian",
        help="Distribution family (now 'gaussian')",
    )
    parser.add_argument("--verbose", action="store_true", help="Print more information")
    parser.add_argument(
        "--n_steps", type=int, default=1000, help="Number of steps for SVI"
    )
    parser.add_argument(
        "--learning_rate",
        "-lr",
        type=float,
        default=0.005,
        help="Learning rate for SVI",
    )
    args = parser.parse_args()
    import warnings

    warnings.filterwarnings("ignore", category=DeprecationWarning)
    main(
        dat_path=args.dat_path,
        locus=args.locus,
        subset_counts=args.subset_counts,
        seed=args.seed,
        max_M=args.max_M,
        holdoff_rate=args.holdoff_rate,
        outDir=args.outDir,
        datName=args.dataset_name,
        no_loop=args.no_loop,
        pick_criterion=args.pick_criterion,
        bfb_family=args.BFB_family,
        verbose=args.verbose,
        n_steps=args.n_steps,
        learning_rate=args.learning_rate,
    )

import pyro
from pyro.distributions.transforms import AffineTransform
from torch.distributions import constraints
import pyro.distributions as dist
import numpy as np
import torch

from mixture_models_pyro_globals import (
    GLOBAL_min_sigma_ln,
    GLOBAL_GET_MIN_SEP,
    GLOBAL_GET_EFFECTIVE_W,
)


# ------------------------------------------------------------------------------
# A custom shifted lognormal that returns -inf when x < loc.
# ------------------------------------------------------------------------------
class ShiftedLogNormal(dist.TransformedDistribution):
    """
    Represents a shifted lognormal distribution.
    If Y ~ LogNormal(log(mu), sigma) and X = Y + loc,
    then the density is
      p(x) = LogNormal(x - loc | log(mu), sigma)
    for x > loc, and 0 otherwise.
    """

    def __init__(self, mu, sigma, loc=0):
        assert loc <= 0, "loc must be <= 0"
        # ln_base = dist.LogNormal(torch.log(mu), sigma, validate_args=False)
        ln_base = dist.LogNormal(mu, sigma, validate_args=False)
        ## shifted
        transform = AffineTransform(loc=loc, scale=1.0)
        super().__init__(ln_base, [transform], validate_args=False)
        self.loc = loc

    def log_prob(self, x):
        mask = x < self.loc
        lp = torch.empty_like(x)
        if (~mask).any():
            lp[~mask] = super().log_prob(x[~mask])
        lp[mask] = -float("inf")
        return lp


# ------------------------------------------------------------------------------
# A custom heterogeneous mixture distribution.
# ------------------------------------------------------------------------------
class HeteroMixture(dist.TorchDistribution):
    """
    A mixture distribution with potentially heterogeneous components.
    The overall density is
      p(x) = sum_i w_i p_i(x)
    where each p_i is a component distribution.
    """

    arg_constraints = {}
    support = constraints.real
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


def make_init_to_custom_value(train_data, M_value, effective_width=None):
    def init_to_custom_value(site):
        name = site["name"]
        if name == "w":
            # return torch.tensor(0.5)
            return torch.tensor(0.5)
        elif name == "log_mu_ln":
            # Initialize to the log of the data mean (adding a small constant for stability)
            return torch.log(torch.tensor(train_data.mean()) + 1e-3)
        elif name == "sigma_ln_orig":
            return torch.tensor(1.0)
        elif name == "loc_ln_free":
            return torch.tensor(0.0)
        elif name == "alpha":
            # Return uniform initial weights
            return torch.ones(M_value) / M_value
        elif name == "mu_components":
            # Initialize each Gaussian mean at the data mean
            # print the value
            print(f"mu_components: {train_data.mean()}")
            return torch.ones(M_value) * train_data.mean()
        elif name == "l_components":
            # Initialize at zero so that sigma_components starts at sigmoid(0)*MAX_GAUSS_SIGMA ~ 0.5*MAX_GAUSS_SIGMA
            return torch.zeros(M_value)
        else:
            # For any other latent variables, fallback to the default initialization
            return site["fn"]()

    return init_to_custom_value


# ------------------------------------------------------------------------------
# The Pyro model.
# ------------------------------------------------------------------------------
def model(data, MAX_GAUSS_SIGMA=None, M=2, effective_width=None):
    """
    Bayesian mixture model with two pieces:
      - A shifted lognormal component with parameters (mu_ln, sigma_ln, loc_ln)
      - A mixture of M Gaussian components with parameters (mu_j, sigma_j)

    The likelihood is:
      p(x) = w * p_shifted_ln(x | mu_ln, sigma_ln, loc_ln)
           + (1-w) * sum_{j=1}^{M} alpha_j * Normal(x | mu_j, sigma_j)
    """
    data = torch.tensor(data, dtype=torch.float)
    N = data.size(0)
    max_val = data.max()

    # --- Hyperparameters for weakly informative priors ---
    a_w, b_w = 1.0, 1.0
    data_mean = data.mean()
    data_std = data.std() + 1e-3
    m0 = torch.log(torch.tensor(2.0))

    # --- Priors for the lognormal component ---
    w = pyro.sample("w", dist.Beta(a_w, b_w))
    log_mu_ln = pyro.sample("log_mu_ln", dist.Normal(m0, 10.0))
    mu_ln = pyro.deterministic("mu_ln", torch.exp(log_mu_ln))

    # Ensure that the data_std is positive
    sigma_ln_orig_sd = np.log(data_std) if data_std > 1 else 1.0
    sigma_ln_orig = pyro.sample("sigma_ln_orig", dist.HalfNormal(sigma_ln_orig_sd))
    sigma_ln = pyro.deterministic("sigma_ln", sigma_ln_orig + GLOBAL_min_sigma_ln)

    loc_ln_free = pyro.sample("loc_ln_free", dist.Normal(0, 10))
    loc_ln = pyro.deterministic("loc_ln", loc_ln_free)
    ln_dist = ShiftedLogNormal(mu_ln, sigma_ln, loc_ln)
    # ln_dist = ShiftedLogNormal(mu_ln, sigma_ln, 0)

    # --- Priors for the Gaussian mixture components ---
    alpha = pyro.sample("alpha", dist.Dirichlet(torch.ones(M)))
    mu0 = pyro.sample("mu0", dist.LogNormal(torch.log(data_mean), 2.0))
    # Define a minimum separation between the means

    # For the remaining M-1 differences, sample positive values.
    delta = pyro.sample("delta", dist.Exponential(1.0).expand([M - 1]).to_event(1))
    # Adjust the differences to guarantee a minimum separation.
    delta_adjusted = delta + GLOBAL_GET_MIN_SEP(effective_width)
    # Construct the ordered Gaussian means.
    mu_components = torch.cat(
        [mu0.unsqueeze(0), mu0.unsqueeze(0) + torch.cumsum(delta_adjusted, dim=0)]
    )

    l_components = pyro.sample(
        "l_components", dist.Normal(0, 1).expand([M]).to_event(1)
    )
    sigma_components = torch.sigmoid(l_components) * MAX_GAUSS_SIGMA

    # --- Softly modulate the weight of the lognormal component ---
    # When max_val is small (<15), soft_factor will be near 0, reducing the lognormal's influence.
    effective_w = GLOBAL_GET_EFFECTIVE_W(w, max_val)

    # --- Build the overall mixture distribution ---
    # comp_weights = torch.cat([w.unsqueeze(0), (1 - w) * alpha])
    comp_weights = torch.cat([effective_w.unsqueeze(0), (1 - effective_w) * alpha])
    gauss_components = [
        dist.Normal(mu_components[i], sigma_components[i]) for i in range(M)
    ]

    components = [ln_dist] + gauss_components
    mix = HeteroMixture(comp_weights, components)

    with pyro.plate("data", N):
        pyro.sample("obs", mix, obs=data)


def model_gaussian(data, MAX_GAUSS_SIGMA=None, M=2, effective_width=None):
    """
    Bayesian mixture model with M Gaussian components that are constrained to be at least

    Sample the first sigma from a .5, .5 beta * the size data.max/6
    If it is small,
        - sample the second sigma as a small one too
        - sample the first and second mu to be anywhere
        - sample alphas from a more concenterated one
    if it is large,
        - sample the second sigma as very small
        - sample the first mu anywhere
        - sample the second mu far apart
        - sample alphas from a more lopsided one

    The likelihood is:

    """
    data = torch.tensor(data, dtype=torch.float)
    N = data.size(0)
    # Sample the first sigma from a .5, .5 beta * the size data.max/6
    sigma0_raw = pyro.sample("sigma0_raw", dist.Beta(0.5, 0.5))
    # Ensure it's between 0 and data.max() / 6
    sigma0 = pyro.deterministic("sigma0", torch.clamp(sigma0_raw * (data.max() / 6), min=0.001, max=data.max() / 6))

    min_delta = effective_width / 2.5 
    # If it is small, i.e., less than 0.1, sample the second sigma as a small one too
    if sigma0 < 0.3:
        # sigma0 is small, so two peaks, just about the same size, 
        # that are apart should explain the data

        # Sample uniform between 0.001 and sigma0
        sigma1_raw = pyro.sample(
            "sigma1_raw", dist.Uniform(0.0, 1.0)
        )
        # Ensure the minimum sigma is at least 0.001
        sigma1 = pyro.deterministic(
            "sigma1", torch.clamp(sigma1_raw * sigma0, min=0.001, max=2.0)
        )
        # Sample the mus whereever
        mu0 = pyro.sample("mu0", dist.Uniform(0, data.max()))
        delta = pyro.sample(
            "delta",
            dist.Uniform(torch.tensor(min_delta), data.max() / 2)
            .expand([M - 1])
            .to_event(1),
        )
        mu1 = pyro.deterministic("mu1", mu0 + delta)
        # Sample alphas from a more concentrated Dirichlet distribution so its more likely evenly
        alphas = pyro.sample("alphas", dist.Dirichlet(5 * torch.ones(M)))
    else:
        # SIGMA0 is large, and should explain the whole distribution
        # So sigma1 should be very smal
        # Make a dummy sigma1_raw here
        
        # ENSURE SIGMA1 is SUPER SMALL i.e., at most .01
        sigma1_raw = pyro.sample(
            "sigma1_raw", dist.Uniform(0.0, 1.0)
        )
        sigma1 = pyro.deterministic(
            "sigma1", torch.clamp(sigma1_raw * .1, min=0.0001, max=1.0)
        )

        # Sample the first mu anywhere
        mu0 = pyro.sample("mu0", dist.Uniform(0, data.max()))
        # Sample the second mu far apart, using a distance variable
        delta = pyro.sample(
            "delta",
            dist.Uniform(torch.tensor(min_delta), data.max() / 2) 
            .expand([M - 1])
            .to_event(1),
        )
        # Adjust the differences to ensure each gap is at least GLOBAL_GET_MIN_SEP(effective_width).
        # Adjust the second mu to be at least delta away from the first mu
        mu1 = pyro.deterministic("mu1", mu0 + delta)
        # Sample alphas from a more lopsided Dirichlet distribution so its more likely to have one component dominate
        # i.e., the first one 1, and the rest 0.5
        alphas = pyro.sample(
            "alphas", dist.Dirichlet(torch.tensor([1.0] + [0.5] * (M - 1)))
        )
    # Create the Gaussian mixture
    # 1. create the mu_component (as torch.cat)
    mu_components = torch.cat(
        [mu0.unsqueeze(0), mu1], dim=0  # shape (1,)  # shape (M-1,)
    )  # result shape (M,)
    # 2. create the sigma_components (as torch.cat)
    sigma_components = torch.cat([sigma0.unsqueeze(0), sigma1.unsqueeze(0)], dim=0)
    # Add to deterministic for easier access later on
    mus = pyro.deterministic("mus_est", mu_components)
    sigmas = pyro.deterministic("sigma_est", sigma_components)
    alphas = pyro.deterministic("alphas_est", alphas)
    # The mixture is built using the mixing weights alpha.
    mix = dist.MixtureSameFamily(
        dist.Categorical(alphas),
        dist.Normal(mu_components, sigma_components),
    )
    with pyro.plate("data", N):
        pyro.sample("obs", mix, obs=data)

def shifted_lognormal_pdf(x, muLN_est, sigmaLN_est, locLn_est):
    """Use the existing shifted lognormal distribution"""
    ln_dist = ShiftedLogNormal(muLN_est, sigmaLN_est, locLn_est)
    # ln_dist = ShiftedLogNormal(muLN_est, sigmaLN_est, 0)
    return ln_dist.log_prob(x).exp()


def gaussian_mixture_pdf(x, alphas_est, mus_est, sigmas_est):
    """Compute the weighted Gaussian mixture density."""
    pdf = torch.zeros_like(x)
    for j in range(len(alphas_est)):
        normal_dist = dist.Normal(mus_est[j], sigmas_est[j])
        pdf += alphas_est[j] * normal_dist.log_prob(x).exp()
    return pdf

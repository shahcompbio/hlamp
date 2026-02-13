import numpy as np
import torch

# Global variables for Pyro Mixture Models
GLOBAL_min_sigma_ln = .5


## Global values
### For gaussian
GLOBAL_gaussian_sigma_min = 0.5 # width --> 1.5
GLOBAL_gaussian_sigma_max = 2.0 # width --> 6.0

### For lognormal

def GLOBAL_GET_MIN_SEP(effective_width):
    """
    Should be based on the effective width of the data.
    The MIN_SEP between Gaussian peaks should not be more than the size of the dataset (effective width)
    The MIN_SEP should be enough to avoid overlap between the two Gaussian peaks.
    """
    min_sep = np.minimum(10.0, 0.5 * effective_width)
    print(f"MIN_SEP: {min_sep} - effective_width: {effective_width}")
    return min_sep


# Global functions
def GLOBAL_GET_MAX_GAUSS_SIGMA(effective_width):
    """
    Get the maximum Gaussian sigma for the mixture model.
    
    Parameters
    ----------
    effective_width : float
        The effective width of the Gaussian.
    
    Returns
    -------
    float
        The maximum Gaussian sigma.
    """
    MAX_GAUSS_SIGMA = np.maximum(0.5 * effective_width / 6.0, 2.0)
    #MAX_GAUSS_SIGMA = np.maximum(0.5 * effective_width / 6.0, 4.0)
    #MAX_GAUSS_SIGMA = np.minimum(MAX_GAUSS_SIGMA, 20.0 / 3.0)
    MAX_GAUSS_SIGMA = np.minimum(MAX_GAUSS_SIGMA, 15.0 / 3.0)
    print(f"MAX_GAUSS_SIGMA: {MAX_GAUSS_SIGMA} - effective_width: {effective_width}")
    return MAX_GAUSS_SIGMA


def GLOBAL_GET_EFFECTIVE_W(w, max_val, max_cn=40.0, scale=2.0):
    """
    Get the effective component weight for the lognormal.
    

    Parameters
    ----------
    w : torch.Tensor
        The weight of the Gaussian.
    max_val : float
        The maximum value in the data.
    
    Returns
    -------
    torch.Tensor
        The effective component weight.
    """
    soft_factor = torch.sigmoid((max_val - max_cn) / scale)
    effective_w = w * soft_factor
    return effective_w


######################
## Math Utils
######################

def get_effective_width(data):
    """ """
    data = np.array(data)
    min_val = np.min(data)
    max_val = np.max(data)
    # Compute a robust range (i.e, remove the outliers)
    q1 = np.percentile(data, 5)
    q99 = np.percentile(data, 95)
    effective_width = q99 - q1
    return effective_width
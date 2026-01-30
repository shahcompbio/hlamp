#!/bin/python
# Run the stochastic agent based simulator for copy number change.

import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d


def plot_population_size(figPath=None, summary=None, final_ndna=None, figsize=(10, 6)):
    plt.clf()
    plt.figure(figsize=figsize)
    plt.plot(summary[:, 0], summary[:, 1], label="Number of Cells")
    plt.title(f"Number of Cells Over Generations {final_ndna.size}")
    plt.xlabel("Generation")
    plt.ylabel("Number of Cells")
    plt.legend()
    if figPath is None:
        plt.show()
    else:
        plt.savefig(figPath, dpi=300, bbox_inches="tight")



def find_tail(data, bins=50):
    mean, std = np.mean(data), np.std(data)
    the_std = np.std(data)
    the_mod = mean + 2 * the_std  # or mean + 3 * std
    # hist, bin_edges = np.histogram(data, bins=bins, density=True)
    # # Compute bin centers
    # bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # # Smooth histogram using a Gaussian filter
    # sigma = 2  # Adjust for more or less smoothing
    # smoothed_hist = gaussian_filter1d(hist, sigma=sigma)
    # Compute the first derivative
    # deriv = np.gradient(smoothed_hist, bin_centers)
    # # Identify the tail threshold (where the derivative significantly drops)
    # tail_index = np.where(deriv > -0.01)[0][-1]  # Last point before flattening out
    # tail_threshold = bin_centers[tail_index]
    #return tail_threshold
    return the_mod


def plot_histogram(figPath=None, final_ndna=None, figsize=(12, 3), show_KDE=False):
    """
    plot a histogram of the number of nDNA
    plot the vector for final_ndna
    plot side by side, the full thing, and then the one with the values above 500

    TODO: Show dashed vertical lines at mean and mod of the distribution
    TODO: Define a way to see the tail of the distribution
    TODO: Highlight the peaks of the distrubiton using the peak-finder module
    """
    nbins = 60
    the_mean = np.mean(final_ndna)
    # define a cutoff for the tail of the distribution (i.e.,?)
    the_mod = find_tail(final_ndna)
    plt.clf()
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    sns.histplot(final_ndna, bins=nbins, kde=show_KDE, ax=axes[0])
    # show the mean and mode of the distribution
    axes[0].axvline(the_mean, color="r", linestyle="--", label="Mean")
    axes[0].axvline(the_mod, color="g", linestyle="--", label="Mode")
    axes[0].set_title("All Values")
    sns.histplot(final_ndna[final_ndna > the_mean], bins=nbins, kde=show_KDE, ax=axes[1])
    axes[1].set_title(f"Values > {the_mean:.2f}")
    sns.histplot(final_ndna[final_ndna > the_mod], bins=nbins, kde=show_KDE, ax=axes[2])
    axes[2].set_title(f"Values > {the_mod:.2f}")
    axes[0].set_xlabel("nDNA Copies")
    axes[1].set_xlabel("nDNA Copies")
    axes[2].set_xlabel("nDNA Copies")
    fig.suptitle("Histogram of nDNA Copies")
    # show the legend just for the first plot
    axes[0].legend()
    plt.tight_layout()
    if figPath is None:
        plt.show()
    else:
        plt.savefig(figPath, dpi=300, bbox_inches="tight")


#def plot_subsample_histogram(figPath=None, final_ndna=None, figsize=(12, 4), n_samples=132, show_KDE=False):
def plot_subsample_histogram(figPath=None, final_ndna=None, figsize=(12, 4), n_samples=1000, show_KDE=False):
    # now randomly sample
    n_samples = min(n_samples, final_ndna.size)
    sampled = np.random.choice(final_ndna, size=n_samples, replace=False)
    plt.figure(figsize=figsize)
    sns.histplot(sampled, bins=40, kde=False)
    # add title
    plt.title(f"Sampled nDNA Copies - {n_samples} Cells")
    plt.xlabel("nDNA Copies")
    plt.ylabel("Counts")
    if figPath is None:
        plt.show()
    else:
        plt.savefig(figPath, dpi=300, bbox_inches="tight")
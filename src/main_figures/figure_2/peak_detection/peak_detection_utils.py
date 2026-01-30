import glob
from glob import glob
import numpy as np
import os
import pandas as pd
from scipy.stats import norm

import numpy as np
import seaborn as sns
import os
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt

## Peack detection specific
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import seaborn as sns
from scipy.signal import peak_widths
import pickle


from mixture_models_pyro_utils import load_data


def setup_dirs(outDir):
    figuresDir = os.path.join(outDir, "figures")
    tablesDir = os.path.join(outDir, "tables")
    dataDir = os.path.join(outDir, "data")
    os.makedirs(dataDir, exist_ok=True)
    os.makedirs(figuresDir, exist_ok=True)
    os.makedirs(tablesDir, exist_ok=True)
    return figuresDir, tablesDir, dataDir



####################################################
# Workhorse peak detection functions
####################################################


def compute_peaks(data, num_bins=None, sigma=1.0, **kwargs):
    if num_bins is None:
        num_bins = np.unique(data).shape[0]
        # minimum of this or the max value in the data
        num_bins = min(num_bins, int(np.max(data)))
    # Padding to capture the left tail
    lo, hi = float(np.min(data)), float(np.max(data))
    bw = (hi - lo) / num_bins
    counts, bin_edges = np.histogram(data, bins=num_bins, range=(lo-bw, hi+bw))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # Gaussian Smoothing using SciPy
    # Adjust sigma for more/less smoothing
    if sigma > 0:
        smoothed_counts_gauss = gaussian_filter1d(counts, sigma)
    else:
        smoothed_counts_gauss = counts
    # Plot the counts histogram
    plt.clf()
    plt.plot(bin_centers, counts, "o--", label="Original Histogram")
    plt.plot(
        bin_centers, smoothed_counts_gauss, label="Smoothed Histogram", linewidth=2
    )
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.title("Histogram and Smoothed Histogram")
    plt.legend()
    plt.show()
    save_path = 'figures/histogram_example.png'
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()
    # Choose one of the smoothed versions for peak detection:
    smoothed_counts = smoothed_counts_gauss
    bin_width = bin_centers[1] - bin_centers[0]
    peaks, properties = find_peaks(smoothed_counts, **kwargs)
    # --- Record the height of each detected peak ---
    peak_heights = smoothed_counts[peaks]
    # --- Compute Base Widths in Index Units ---
    # Use rel_height=1.0 to measure the width of each peak.
    results_base = peak_widths(smoothed_counts, peaks, rel_height=1.0)
    widths_indices = results_base[0]  # widths in index units
    left_ips = results_base[2]  # left interpolated positions (indices)
    right_ips = results_base[3]  # right interpolated positions (indices)
    # --- Convert Widths from Index Units to x-Axis Units ---
    widths_x = widths_indices * bin_width
    results_dict = {
        "bin_centers": bin_centers,
        "counts": counts,
        "smoothed_counts": smoothed_counts,
        "peaks": peaks,
        "properties": properties,
        "peak_heights": peak_heights,
        "widths_x": widths_x,
        "left_ips": left_ips,
        "right_ips": right_ips,
        "bin_width": bin_width,
        "widths_indices": widths_indices,
    }
    return results_dict


def process_peak_results(
    bin_centers, peaks, left_ips, right_ips, bin_width, widths_indices, widths_x
):
    """
    Compute the left_base_x, and right_base_x for each peak
    """
    left_base_x = bin_centers[0] + left_ips[i] * bin_width
    right_base_x = bin_centers[0] + right_ips[i] * bin_width
    return left_base_x, right_base_x


def peak_driver(data):
    """
    Given the data, i.e., a list of CN values, comptue the left and right base positions of the peaks
    """
    results = compute_peaks(data)
    bin_centers = results["bin_centers"]
    peaks = results["peaks"]
    left_ips = results["left_ips"]
    right_ips = results["right_ips"]
    bin_width = results["bin_width"]
    left_base_positions = []
    right_base_positions = []
    for i in range(len(peaks)):
        left_base_x = bin_centers[0] + left_ips[i] * bin_width
        right_base_x = bin_centers[0] + right_ips[i] * bin_width
        # Convert to flat python values
        left_base_x = float(left_base_x)
        right_base_x = float(right_base_x)
        left_base_positions.append(left_base_x)
        right_base_positions.append(right_base_x)
    return left_base_positions, right_base_positions


def print_resutls(
    bin_centers, peaks, left_ips, right_ips, bin_width, widths_indices, widths_x
):
    # --- Print the Results ---
    for i, peak in enumerate(peaks):
        print(f"Peak at x = {bin_centers[peak]:.2f}:")
        print(f"  Base width (index units): {widths_indices[i]:.2f}")
        print(f"  Base width (x-axis units): {widths_x[i]:.2f}")
        # For approximate positions of the left and right bases in x-units:
        left_base_x = bin_centers[0] + left_ips[i] * bin_width
        right_base_x = bin_centers[0] + right_ips[i] * bin_width
        print(
            f"  Left base at x = {left_base_x:.2f}, Right base at x = {right_base_x:.2f}\n"
        )


def plot_result(
    bin_centers,
    counts,
    smoothed_counts,
    peaks,
    left_ips,
    right_ips,
    bin_width,
    figsize=(6, 6),
    fpath=None,
    titleStr=None,
    data=None,
):
    """
    Plot side by side with a histogram
    """
    # --- Plotting ---
    plt.clf()
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    # On the right: histogram and smoothed histogram
    ax = axes[0]
    ax.plot(bin_centers, counts, "o--", label="Original Histogram")
    ax.plot(bin_centers, smoothed_counts, label="Smoothed Histogram", linewidth=2)
    ax.plot(
        bin_centers[peaks],
        smoothed_counts[peaks],
        "rx",
        markersize=10,
        label="Detected Peaks",
    )
    # Optionally, plot horizontal lines indicating the base widths
    for i in range(len(peaks)):
        left_base_x = bin_centers[0] + left_ips[i] * bin_width
        right_base_x = bin_centers[0] + right_ips[i] * bin_width
        # Determine the y level at which the base is measured. Typically, this is the height of the peak
        # minus its prominence down to the baseline. Here, for simplicity, we use the minimal of the
        # values at the interpolated left/right positions (if meaningful).
        base_y = min(
            np.interp(left_ips[i], np.arange(len(smoothed_counts)), smoothed_counts),
            np.interp(right_ips[i], np.arange(len(smoothed_counts)), smoothed_counts),
        )
        ax.hlines(
            y=base_y,
            xmin=left_base_x,
            xmax=right_base_x,
            color="green",
            linestyle="--",
            label="Base Width" if i == 0 else "",
        )
    # set the x and y labels
    ax.set_xlabel("Value")
    ax.set_ylabel("Frequency")
    # plt the left: histogram only
    ax = axes[1]
    # use sns
    sns.histplot(data, bins=len(counts), kde=False, alpha=0.5, ax=ax)
    # show the peaks on the true data (i.e., x: data value at the peak, y: the peak value) - use sns
    sns.scatterplot(
        x=bin_centers[peaks],
        y=smoothed_counts[peaks],
        ax=ax,
        color="red",
        s=100,
        label="Detected Peaks",
    )
    if titleStr is not None:
        plt.title(titleStr)
    else:
        plt.title("Histogram, Peak Detection, and Base Width Calculation")
    plt.legend()
    if fpath is None:
        plt.show()
    else:
        plt.savefig(fpath, dpi=300, bbox_inches="tight")


def handle_data(
    dat,
    outDir,
    verbose=False,
    fname=None,
    sigma=1,
    num_bins=None,
    do_plot=True,
    **kwargs,
):
    results = compute_peaks(dat, sigma=sigma, num_bins=num_bins, **kwargs)
    # create outDir
    os.makedirs(outDir, exist_ok=True)
    bin_centers = results["bin_centers"]
    counts = results["counts"]
    smoothed_counts = results["smoothed_counts"]
    peaks = results["peaks"]
    left_ips = results["left_ips"]
    right_ips = results["right_ips"]
    bin_width = results["bin_width"]
    widths_x = results["widths_x"]
    widths_indices = results["widths_indices"]
    left_base_positions = []
    right_base_positions = []
    peak_heights_list = []
    for i in range(len(peaks)):
        left_base_x = bin_centers[0] + left_ips[i] * bin_width
        right_base_x = bin_centers[0] + right_ips[i] * bin_width
        # Convert to flat python values
        left_base_x = float(left_base_x)
        right_base_x = float(right_base_x)
        left_base_positions.append(left_base_x)
        right_base_positions.append(right_base_x)
        # collect heights
        peak_heights_list.append(float(results["peak_heights"][i]))
    if verbose:
        print_resutls(
            bin_centers, peaks, left_ips, right_ips, bin_width, widths_indices, widths_x
        )
    fpath = os.path.join(outDir, fname)
    titleStr = fname.split(".")[0]
    # add the max of width_x to the fname
    print(f"width_x len is {len(widths_x)}")
    if len(widths_x) > 0:
        max_width = np.max(widths_x)
        # replace the .png with _max_width.png
        print(f"Max width is {max_width}")
        fpath = fpath.replace(".png", f"_{max_width:.2f}.png")
    if do_plot:
        plot_result(
            bin_centers,
            counts,
            smoothed_counts,
            peaks,
            left_ips,
            right_ips,
            bin_width,
            figsize=(10, 6),
            fpath=fpath,
            titleStr=titleStr,
            data=dat,
        )
    return left_base_positions, right_base_positions, peak_heights_list, results


def handle_item(
    fpath=None,
    locus=None,
    tumor_cells=None,
    figuresDir=None,
    dataDir=None,
    drop_normals=True,
    minimal_normal=5,
    num_bins=100,
    **kwargs,
):
    """
    Computes peaks for a given fpath and locus
    
    Parameters
    ----------
    num_bins : int
        The more, the more peaks to be found. Default is 100 which skips the left tail peaks.
        
    """
    assert figuresDir is not None, "Please provide figuresDir"
    assert dataDir is not None, "Please provide dataDir"
    # Load the data
    dataset_name = os.path.basename(fpath).split(".")[0]
    vals, _, _, remaining_cells = load_data(
        dat_path=fpath,
        locus=locus,
        do_log=False,
        do_prune_left=False,
        keep_cells=tumor_cells,
        drop_normals=drop_normals,
        minimal_normal=minimal_normal,
    )
    vals_df = pd.DataFrame({"cell_id": remaining_cells, "value": vals})
    # Print the minium and maximum values
    print(f"Dataset: {dataset_name}, Locus: {locus}, Min value: {np.min(vals)}, Max value: {np.max(vals)}, N cells: {vals.shape[0]}")
    n_cells = vals.shape[0]
    fname = f"{dataset_name}_{locus}_peak_detection.png"
    left_base_positions, right_base_positions, peak_heights_list, orig_results = (
        handle_data(
            vals,
            num_bins=100,
            sigma=1.5,
            outDir=figuresDir,
            fname=fname,
            height=1.5,
            distance=1,
            width=(2, None),
            prominence=1.5,
            # pass on the options
            **kwargs,
        )
    )
    # Save the results as a pickle files
    results_dict = {
        "left_base_positions": left_base_positions,
        "right_base_positions": right_base_positions,
        "peak_heights_list": peak_heights_list,
        "bin_centers": vals,
        "dataset_name": dataset_name,
        "locus": locus,
        "n_cells": n_cells,
        "orig_results": orig_results,
        "vals_df": vals_df,
    }
    pickle_path = os.path.join(dataDir, f"{dataset_name}_{locus}_peak_detection.pkl")
    with open(pickle_path, "wb") as f:
        pickle.dump(results_dict, f)
    n_peaks = len(left_base_positions)
    return dataset_name, locus, n_peaks, pickle_path


def get_paths(individual=None, ignore_qc=False, qc_files=None, oncogene_files=None):
    assert oncogene_files is not None, "Please provide oncogene_files"
    fpath = [f for f in oncogene_files if individual in f]
    assert (
        len(fpath) == 1
    ), f"Expected one file for individual {individual}, found {len(fpath)}"
    fpath = fpath[0]
    if ignore_qc:
        qc_path = None
    else:
        qc_path = [f for f in qc_files if individual in f]
        assert (
            len(qc_path) == 1
        ), f"Expected one QC file for individual {individual}, found {len(qc_path)}"
        qc_path = qc_path[0]
    return fpath, qc_path


def get_tumor_cells(qc_path=None):
    qc_df = pd.read_csv(qc_path)
    # Only keep cell_ids where the is_tumor_cells == 1
    tumor_cells = qc_df[qc_df["is_tumor_cell"] == "yes"]["cell_id"].tolist()
    return tumor_cells

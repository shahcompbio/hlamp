import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.ndimage import gaussian_filter1d
from scipy.signal import peak_widths


def compute_peaks(data, num_bins=None, sigma=1.0, **kwargs):
    """
    TODO: this function will throw an error if there is only one peak found,
    e.g., when the max value in data is very small
    """
    if num_bins is None:
        num_bins = np.unique(data).shape[0]
        # minimum of this or the max value in the data
        num_bins = min(num_bins, int(np.max(data)))
    if num_bins < 2:
        print(f"Found {num_bins}, will return None")
        return None
    counts, bin_edges = np.histogram(data, bins=num_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # Gaussian Smoothing using SciPy
    # Adjust sigma for more/less smoothing
    if sigma > 0:
        smoothed_counts_gauss = gaussian_filter1d(counts, sigma)
    else:
        smoothed_counts_gauss = counts
    # Choose one of the smoothed versions for peak detection:
    smoothed_counts = smoothed_counts_gauss  # or use smoothed_counts_ma
    bin_width = bin_centers[1] - bin_centers[0]
    peaks, properties = find_peaks(smoothed_counts, **kwargs)
    # --- Compute Base Widths in Index Units ---
    # Use rel_height=1.0 to measure the width at the base of each peak.
    results_base = peak_widths(smoothed_counts, peaks, rel_height=1.0)
    widths_indices = results_base[0]   # widths in index units
    left_ips = results_base[2]         # left interpolated positions (indices)
    right_ips = results_base[3]        # right interpolated positions (indices)
    # --- Convert Widths from Index Units to x-Axis Units ---
    widths_x = widths_indices * bin_width
    results_dict = {
        'bin_centers': bin_centers,
        'counts': counts,
        'smoothed_counts': smoothed_counts,
        'peaks': peaks,
        'properties': properties,
        'widths_x': widths_x,
        'left_ips': left_ips,
        'right_ips': right_ips,
        'bin_width': bin_width,
        'widths_indices': widths_indices
    }
    return results_dict


def process_peak_results(bin_centers, peaks, left_ips, right_ips, bin_width, widths_indices, widths_x):
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
    if results is None:
        return [], []
    bin_centers = results['bin_centers']
    peaks = results['peaks']
    left_ips = results['left_ips']
    right_ips = results['right_ips']
    bin_width = results['bin_width']
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


def load_data(
    dat_path,
    locus=None,
    do_log=False,
    do_inject_ref=False,
    do_inject_tail=False,
    do_prune_left=True,
    keep_cells=None,
    drop_normals=True,
    minimal_normal=5,
):
    """
    Load HMMCopy processed data.
    Assumes the data is in the format of:
    - columns: the loci, listed by name
    - rows: the cells

    Args:
    -----
        do_inject_ref: bool, whether to add artificial cells with copy number 2.0
        do_inject_tail: bool, whether to inject the tail values above the 95% percentile
        do_cutoff: bool, whether to cut off values below 10
        keep_cells: list of cells to keep, if None, keep all cells
        drop_normals: remove all entries with copy number 5 or below
    """
    sep = "," if (dat_path.endswith(".csv") or dat_path.endswith(".csv.gz")) else "\t"
    dat = pd.read_csv(dat_path, sep=sep)
    # filter for keep_cells if provided
    if keep_cells is not None:
        print(f"Filtering for {len(keep_cells)} cells. Original data has {dat.shape[0]} cells")
        dat = dat[dat['cell_id'].isin(keep_cells)].copy()
        if dat.shape[0] == 0:
            print("No cells found after filtering, returning empty array")
            return np.array([]), None, np.array([])
    # print(dat.head())
    assert locus is not None, "Please provide a locus to extract"
    assert locus in dat.columns, "Locus not found in data"
    vals = dat[locus].values
    orig_shape = vals.shape
    # Remove nans
    vals = vals[~np.isnan(vals)]
    shape_sans_nans = vals.shape
    # Remove negative values
    vals = vals[vals >= 0]
    shape_sans_negs = vals.shape
    all_vals = vals.copy()
    print(
        f"Loaded data with shape {orig_shape}, removed nans to get {shape_sans_nans}, removed negatives to get {shape_sans_negs}"
    )
    # inject two observations at 2.0: to anchor the references
    if do_inject_ref:
        n_ref_inject = 5
        extra_refs = np.array([2.0] * n_ref_inject)
        vals = np.concatenate([vals, extra_refs])
        print(f"Injecting {n_ref_inject} extra references at 2.0")
        vals = np.concatenate([vals, np.array([2.0, 2.0])])
    tail_vals = None
    if do_inject_tail:
        # Inject the 95% percentile values or more to the tail
        ninety_five_percentile = np.percentile(vals, 95)
        # Extract all values above the 95% percentile
        tail_vals = vals[vals >= ninety_five_percentile]
        print(f"Injecting tail values above 95% percentile: {tail_vals}")
    if do_prune_left:
        # Remove values that are smaller than the mean frequncy and have copy number less than 10
        left_base_positions, right_base_positions = peak_driver(vals)
        # if left_base_positions is not empty, and its value is less than 10, remove all values to the left of it
        if len(left_base_positions) > 0:
            print(f"Left base positions: {np.round(left_base_positions, 2)}")
            for left_base_position in left_base_positions:
                if left_base_position < 10:
                    # Remove all the values to the left of it
                    print(f"Removing all values to the left of {left_base_position}")
                    vals = vals[vals >= left_base_position].copy()
                else:
                    # Otherwise, leave
                    break
    if drop_normals:
        print(f"Dropping all values with copy number {minimal_normal} or below")
        vals = vals[vals > minimal_normal].copy()
        print(f"Remaining cells: {len(vals)}")
    if do_log:
        print("Applying log transform")
        vals = np.log(vals + 1e-10)
        all_vals = np.log(all_vals + 1e-10)
        len_before = len(vals)
        # Remove all negative values
        vals = vals[vals >= 0.0].copy()
        len_after = len(vals)
        print(f"Cells remaining {len_after}/{len_before}")
        if tail_vals is not None:
            # Apply log transform to tail values as well
            tail_vals = np.log(tail_vals + 1e-10)
            tail_vals = tail_vals[tail_vals >= 0.0].copy()
    return vals, tail_vals, all_vals

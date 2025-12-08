import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager
from matplotlib import rcParams
from sklearn.preprocessing import StandardScaler
from mixture_models_pyro_utils import load_data
import yaml
from scipy.stats import norm

rcParams["font.family"] = "Arial"


# Compute a weighted sum
def compute_weighted_sum(vals):
    """
    What is the effective copy number?
    Effective copy number (ECN):
        - w = sum(val * weight for val, weight in vals)

    Args:
    -----
    vals: a list of copy numbers
    """
    # 1. compute the weight (i.e., freq) of each copy number bin
    n_bins = 100
    # bin the data
    counts, bin_edges = np.histogram(vals, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    # 2. normalize counts to compute weights
    weights = counts / counts.sum()
    # 3. compute the weighted sum
    weighted_sum = sum(center * weight for center, weight in zip(bin_centers, weights))
    return weighted_sum


def make_01(x):
    """Normalize a numpy array to [0, 1]"""
    return (x - x.min()) / (x.max() - x.min())

def compute_metrics(res, mean_coeff=0.5, alpha_coeff=.5):
    x = res["sigma_alpha_max"].values
    y = res["max_alpha"].values
    # standardize them
    # scaler = StandardScaler()
    # x_norm = scaler.fit_transform(x.reshape(-1, 1)).flatten()
    # y_norm = scaler.fit_transform(y.reshape(-1, 1)).flatten()
    # # normalize to [0,1]
    x_norm = make_01(x)
    y_norm = make_01(y)
    # compute d_parallel and d_perp in one go:
    d_par = (x_norm + y_norm) / np.sqrt(2)
    d_perp = np.abs(y_norm - x_norm) / np.sqrt(2)
    res["d_parallel"] = np.max(d_par) - d_par
    res["d_perpendicular"] = d_perp
    res["w2"] = res["d_parallel"] * alpha_coeff + res["d_perpendicular"] * (1 - alpha_coeff)
    # Compute w3 as the max(w2) - w2
    res["w3"] = res["w2"].max() - res["w2"]
    # Normalize between 0 and 1
    res["w3"] = (res["w3"] - res["w3"].min()) / (res["w3"].max() - res["w3"].min())
    # sort by w3
    res.sort_values(by="w3", ascending=False, inplace=True)
    # Define w4 as a metric that combines w3 with mu_alpha_max/100
    # z-score mu_alpha_max using standard scaler
    scaler = StandardScaler()
    #res["mu_alpha_max_z"] = scaler.fit_transform(res[["mu_alpha_max"]])
    res["mu_alpha_max_z"] = res["mu_alpha_max"].values
    # between 0 and 1
    res["mu_alpha_max_z"] = make_01(res["mu_alpha_max_z"])
    #res["w4"] = res["w3"] + mean_coeff * res["mu_alpha_max_z"] * res["max_alpha"]
    res["w4"] = res["w3"] + mean_coeff * res["mu_alpha_max_z"] * res["max_alpha"]
    # Between 0 and 1
    res["w4"] = make_01(res["w4"])
    res.sort_values(by="w4", ascending=False, inplace=True)
    return res



#####################################################
## Loop utils
#####################################################

def gather_best_model_results(paths):
    """
    Paths to best_model.yaml files
    Usually under **/**/best_model.yaml
    """
    # only keep unique
    paths = list(set(paths))
    res = None
    for pp in paths:
        try:
            # load the yaml file
            with open(pp, 'r') as file:
                config = yaml.safe_load(file)
            mm_path = pp.replace('best_model.yaml', 'configs.yaml')
            with open(mm_path, 'r') as file:
                config_mm = yaml.safe_load(file)
            # get the best model
            system_name = pp.split('/')[-3]
            locus = pp.split('/')[-2]
            # get the alphas
            alphas = config['alphas_est']
            # get the max alpha
            max_alpha = max(alphas)
            # get the min alpha
            min_alpha = min(alphas)
            # compute the diff 
            diff = max_alpha - min_alpha
            # also find sigma_alpha_max and sigma_alpha_min
            sigmas = config['sigmas_est']
            sigma_alpha_max = sigmas[alphas.index(max_alpha)]
            sigma_alpha_min = sigmas[alphas.index(min_alpha)]
            # for mu
            mu = config['mus_est']
            mu_alpha_max = mu[alphas.index(max_alpha)]
            mu_alpha_min = mu[alphas.index(min_alpha)]
            tmp = {}
            tmp['system_name'] = system_name
            tmp['locus'] = locus
            tmp['max_alpha'] = max_alpha
            tmp['min_alpha'] = min_alpha
            tmp['diff'] = diff
            tmp['sigma_alpha_max'] = sigma_alpha_max
            tmp['sigma_alpha_min'] = sigma_alpha_min
            tmp['mu_alpha_max'] = mu_alpha_max
            tmp['mu_alpha_min'] = mu_alpha_min
            # Add the mass as well
            # {'x_mode': 3.0400664806365967, 'density_at_mode': 0.7118465304374695, 'mass_in_window': 0.6161591410636902}
            tmp['x_mode'] = config.get('x_mode', None)
            tmp['density_at_mode'] = config.get('density_at_mode', None)
            tmp['mass_in_window'] = config.get('mass_in_window_0.5', None)
            # Add the config_mm info as well
            # load the data and compute its sd
            vals, _, _ = load_data(config_mm['dat_path'], locus=config_mm['locus'], do_log=True)
            # Fit a Gaussian and extract the sd
            mu, sd = norm.fit(vals)
            ecn = compute_weighted_sum(vals)
            tmp['ecn'] = ecn
            tmp['normal_sd'] = sd
            # Compute the sd
            tmp['max_cn'] = vals.max()
            tmp['sd'] = vals.std()
            # just keep .2 decimals
            for key in tmp.keys():
                if isinstance(tmp[key], float):
                    tmp[key] = round(tmp[key], 3)
            tmp['effective_width'] = config_mm['effective_width']
            tmp['n_cells'] = config_mm['n_cells']
            tmp_df = pd.DataFrame(tmp, index=[0])
            if res is None:
                res = tmp_df
            else:
                res = pd.concat([res, tmp_df], axis=0)
        except Exception as e:
            print(f"Error processing {pp}: {e}")
            continue
    return res


#####################################################
## Plotting
#####################################################

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text


def plot_scatter(res, plot_col=None, min_n_cells=20, n_top=30, fpath=None, dashed_line=None):
    assert plot_col is not None, "plot_col must be specified"
    assert plot_col in res.columns, f"plot_col {plot_col} not found in res columns"
    # 1. Filter your data
    if 'id' not in res.columns:
        res['id'] = res['system_name'].astype(str) + '_' + res['locus'].astype(str)
    res_sub = res[res['n_cells'] > min_n_cells].copy()
    # 2. Create a numeric 'id_code' in the order each ID first appears
    res_sub['id_code'], uniques = pd.factorize(res_sub['id'])
    # 3. Plot
    fig, ax = plt.subplots(figsize=(50, 10))
    sns.scatterplot(
        data=res_sub,
        x="id_code",
        y=plot_col,
        ax=ax,
        size=10,
        linewidth=0,
        alpha=0.7,
    )
    if dashed_line is not None:
        # Add a horizontal dashed line at the specified value
        ax.axhline(y=dashed_line, color='red', linestyle='--', label=f'{dashed_line}')
    # Titles and labels
    ax.set_title(f"{plot_col} vs ID (n_cells > {min_n_cells}) (#items= {res_sub.shape[0]})")
    ax.set_xlabel("ID")
    ax.set_ylabel(plot_col)
    # 4. Annotate the top 20 by plot_col
    top20 = res_sub.nlargest(n_top, plot_col)
    texts = []
    for _, row in top20.iterrows():
        texts.append(
            ax.text(
                row['id_code'],
                row[plot_col],
                row['id'],
                fontsize=8
            )
        )
    # # 5. Adjust text with arrows
    _ = adjust_text(
        texts,
        arrowprops=dict(arrowstyle='->', color='black', lw=0.5),
        #force_points=0.1,
        #expand_points=(1.2, 1.2), 
        # ensure the text is on the right side of the arrow
    )
    # 6. (Optional) Restore x‐axis tick labels in original ID order
    _ = ax.set_xticks(range(len(uniques)))
    _ = ax.set_xticklabels(uniques, rotation=90, fontsize=6)
    ax.tick_params(axis="x", bottom=False, top=False)
    # 7. Legend, layout, and save
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    # set limits to be between 0 and 1.2
    ax.set_ylim(0, 1.2)
    # rmeove the grid
    ax.grid(False)
    fig.savefig(
        fpath,
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(fig)


def plot_metrics(res=None, fpath=None):
    """
    Plot these on a scatter plot (X: sigma_alpha_max, Y: max_alpha, color: [w2, d_parapplel, d_perpendicular],
    """
    plt.clf()
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    # For each meteric, plot the scatter plot
    metrics = ["w2", "d_parallel", "d_perpendicular"]
    for i, metric in enumerate(metrics):
        ax = axes[i]
        sns.scatterplot(
            data=res,
            x="sigma_alpha_max",
            y="max_alpha",
            hue=metric,
            palette="viridis",
            ax=ax,
            linewidth=0,
            alpha=0.7,
        )
        ax.set_title(f"{metric}")
        ax.set_xlabel("Sigma Alpha Max")
        ax.set_ylabel("Max Alpha")
        # ax.legend(title=metric, loc='upper left')
        # rmeove the legend
        ax.legend([], [], frameon=False)
        ax.grid(False)
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    # Save the figure
    fig.savefig(
        fpath,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


# plot on a scatter plot (x: effective_width, y: w_est, color: group)
def plot_embedding(res, min_cells=15, fpath=None):
    df = res.copy()
    df = df[df["n_cells"] > min_cells]
    df["n_cells"] = df["n_cells"]
    x_col_names = ["sigma_alpha_max"]
    i = 0
    df["group"].value_counts()
    df[df["group"] == "ecDNA"]
    # set nan to 'NA'
    df["group"] = df["group"].fillna("no")
    color_dict = {"ecDNA": "orange", "no": "blue"}
    # set dpi to 300
    plt.rcParams["figure.dpi"] = 300
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    # plot the x_col_names on the x axis
    sns.scatterplot(
        data=df,
        x=x_col_names[i],
        y="max_alpha",
        hue="group",
        palette=color_dict,
        size="n_cells",
        ax=ax,
    )
    # set the title
    ax.set_title(f"#cells > {min_cells}  (#items= {df.shape[0]})")
    # set the x and y labels
    ax.set_xlabel("Width of largest peak (sigma_alpha_max)")
    ax.set_ylabel("max_alpha")
    # set the x and y limits
    # set the legend
    ax.legend(title="Group", loc="upper left")
    # set the grid
    ax.grid(False)
    # put the legend outside the plot
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.show()
    fig.savefig(fpath, dpi=300, bbox_inches="tight")
    plt.close(fig)

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import scgenome
from scgenome.tools.ranges import dataframe_to_pyranges, pyranges_to_dataframe
import pandas as pd
import warnings
import anndata as ad
import matplotlib.colors as mcolors
from collections import defaultdict

# global constant gene boundaries for genes of interest
gene_bounds = defaultdict(lambda:[])
gene_bounds['1'] = [('MDM4', [204485511, 204542871])]
gene_bounds['2'] = [('MYCN', [16080686, 16087129])]
gene_bounds['5'] = [('TERT', [1253282, 1295183])]
gene_bounds['7'] = [('EGFR', [55086710, 55279321]),
                   ('ZNF713', [55955149, 56009918])]
gene_bounds['8'] = [('MYC', [128747680, 128753674])]
gene_bounds['11'] = [('CCND1', [69455924, 69469242])]
gene_bounds['12'] = [('CDK4', [58141510, 58146093]),
                    ('MDM2', [69201952, 69244466])]

# post-ILP processing
def compute_residuals(data, loadings, amplicons, absolute=True, return_values=False):
    n_amplicons = amplicons.shape[0]
    values = np.array(
                [[sum(loadings[i, p] * amplicons[p, j] for p in range(n_amplicons))
                for i in range(data.shape[0])] for j in range(data.shape[1])]).T
    residuals = data - values
    
    if absolute:
        if return_values:
            return values, np.abs(residuals) 
        else:
            return np.abs(residuals)
    else:
        if return_values:
            return values, residuals
        else:
            return residuals

def compose_results(hq_adata, amplicon_adata, data, result):
    # same rows as hq_adata, same columns as amplicon_adata
    result_adata = amplicon_adata[hq_adata.obs.index].copy()

    amplicon_starts = [(np.min(np.where(result['amplicons'][i] > 0)[0]) if np.any(result['amplicons'][i] > 0) else np.inf) for i in range(result['amplicons'].shape[0])]
    amplicon_order = np.argsort(amplicon_starts)
    result['amplicons'] = result['amplicons'][amplicon_order]
    result['loadings'] = result['loadings'][:, amplicon_order]
    
    # remove amplicons with empty loadings
    empty_loadings = np.where(np.sum(result['loadings'], axis=0) == 0)[0]
    result['amplicons'] = np.delete(result['amplicons'], empty_loadings, axis=0)
    result['loadings'] = np.delete(result['loadings'], empty_loadings, axis=1)
    
    result_adata.layers['data'] = data
    values, residuals = compute_residuals(data, result['loadings'], result['amplicons'], absolute=False, return_values=True)
    result_adata.layers['residual'] = residuals
    result_adata.layers['fitted_values'] = values
    result_adata.obs['total_abs_resid'] = np.sum(np.abs(result_adata.layers['residual']),  axis=1)
    result_adata.obs['total_magnitude'] = np.sum(np.abs(result_adata.layers['data']),  axis=1)
    result_adata.obs = result_adata.obs.merge(hq_adata.obs[hq_adata.obs.columns.difference(result_adata.obs.columns)], 
                        left_index=True, right_index=True, how='left')

    result_adata.obsm['loadings'] = result['loadings']
    result_adata.varm['amplicons'] = result['amplicons'].T

    # cluster cells based on loadings
    best_bic = np.inf
    best_clustering = None
    for k in range(2, 10):
        clustering, bic = scgenome.tools.cluster._gmm_diag_bic(result_adata.obsm['loadings'], k)
        if bic < best_bic:
            best_clustering = clustering
            best_bic = bic
    
    result_adata.obs['cluster_id'] = best_clustering
    result_adata.obs['clustering_order'] = np.argsort(best_clustering)
    return result_adata

def plot_amplicons(result_adata, realval_amplicons=True, ax=None, legend=True, include_shading=True,
marker='.', s=None, offset_factor=20):
    if ax is None:
        ax = plt.gca()
        
    amplicons = result_adata.varm['amplicons'].T
    if amplicons.shape[0] > 10:
        cmap = plt.get_cmap('tab20')
    else:
        cmap = plt.get_cmap('tab10')
    if realval_amplicons:
        
        for i in range(amplicons.shape[0]):
            if offset_factor > 0:
                ax.plot(amplicons[i] + i/offset_factor, c=cmap(i/10), label=f'amplicon {i}')
            else:
                ax.plot(amplicons[i], c=cmap(i/10), label=f'amplicon {i}')
        ax.set_xlabel('bin')
        ax.set_ylabel('Number of copies')
    else:
        for i in range(amplicons.shape[0]):
            nz = np.where(amplicons[i] > 0)[0]
            ax.scatter(nz, [i] * len(nz), color=cmap(i/10), label=f'amplicon {i}', marker=marker, s=s)
        ax.set_xlabel('10kb bin')
        ax.set_ylabel('Amplicon index')

    _, ymax = ax.get_ylim()
    boundaries = []
    for i, (amplicon_id, avar) in enumerate(result_adata.var.groupby('amplicon_id', observed=True)):
        start = avar.iindex.min() - 0.5
        end = avar.iindex.max() + 0.5
        if include_shading:
            ax.add_patch(Rectangle((start, 0), 
                            width=(end-start), 
                            height=ymax, zorder=-1, label=amplicon_id, facecolor = plt.get_cmap('tab20')(i), alpha=0.5))
            for discontinuity in np.where(np.diff(avar.iindex) != 1)[0]:
                ax.axvline(x=discontinuity + start + 1, color='grey', linestyle='--')
        boundaries.append(end)
        
    if legend:
        ax.legend()
        sns.move_legend(ax, loc='upper left', bbox_to_anchor=(1,1))
    return boundaries

def plot_residuals(result_adata, residual_bounds=[-20, 20], data_vmax=1000, loading_vmax=None, realval_amplicons=True, relative_resid=False):
    plt.figure(figsize=(20,3), dpi=300)
    boundaries = plot_amplicons(result_adata, realval_amplicons=realval_amplicons) 
        
    plt.figure(figsize=(12,6), dpi=300)
    plt.subplot(1, 3, 1)
    sns.heatmap(result_adata.layers['data'][result_adata.obs['clustering_order']], vmax=data_vmax, yticklabels=[])
    for b in boundaries:
        plt.gca().axvline(x=b, color='cyan', linewidth=0.6)
    plt.title("data")
    
    if loading_vmax is None:
        loading_vmax = np.max(result_adata.obsm['loadings'])
    plt.subplot(1, 3, 2)
    sns.heatmap(result_adata.obsm['loadings'][result_adata.obs['clustering_order']], yticklabels=[], vmax=loading_vmax)
    plt.title('loadings')
    
    plt.subplot(1, 3, 3)
    if relative_resid:
        plotmat = result_adata.layers['residual'][result_adata.obs['clustering_order']] / np.maximum(result_adata.layers['data'], 1)
        plt.title('relative error')
    else:
        plotmat = result_adata.layers['residual'][result_adata.obs['clustering_order']]
        plt.title('residuals')
    sns.heatmap(plotmat, yticklabels=[], cmap='coolwarm', vmin=residual_bounds[0], 
                vmax=residual_bounds[1])
    for b in boundaries:
        plt.gca().axvline(x=b, color='cyan', linewidth=0.6)
    
    
    plt.tight_layout()


def add_cn_bin(adata, cn_bins, id_field='snv_id'):
    assert "bin" in cn_bins

    # Convert from HMMCopy 1-based end included to 1-based end excluded pyranges
    cn_ranges = dataframe_to_pyranges(cn_bins)
    cn_ranges = cn_ranges.assign("End", lambda df: df["End"] + 1)

    snv_ranges = dataframe_to_pyranges(
        adata.var.reset_index().assign(
            chr=lambda df: df["chr"],
            start=lambda df: df["pos"],
            end=lambda df: df["pos"],
        )
    )
    snv_ranges = snv_ranges.assign(
        "End", lambda df: df["End"] + 1
    ) 

    i1 = snv_ranges.intersect(cn_ranges).as_df()
    i2 = cn_ranges.intersect(snv_ranges).as_df()
    intersection = i1.merge(i2)

    adata.var["cn_bin"] = intersection.set_index(id_field)["bin"]

    #assert not adata.var["cn_bin"].isnull().any()

    return adata


def plot_cell_fit(result_adata, cell_id, data_layer_name='data', value_layer_name='fitted_values', residual_layer_name='residual'):
    idx = result_adata.obs.index.get_loc(cell_id)

    plt.plot(result_adata.layers['data'][idx], label='data', color='black')
    plt.plot(result_adata.layers['fitted_values'][idx], label='fitted values', color='grey')
    
    ymin, ymax = plt.ylim()
    boundaries = []
    for i, (amplicon_id, avar) in enumerate(result_adata.var.groupby('amplicon_id', observed=True)):
        start = avar.iindex.min() - 0.5
        end = avar.iindex.max() + 0.5
        plt.gca().add_patch(Rectangle((start, ymin), 
                           width=(end-start), 
                           height=(ymax-ymin), zorder=-1, label=amplicon_id, facecolor = plt.get_cmap('tab20')(i), alpha=0.5))
        for discontinuity in np.where(np.diff(avar.iindex) != 1)[0]:
            plt.gca().axvline(x=discontinuity + start + 1, color='grey', linestyle='--')
        boundaries.append(end)
        
    plt.legend()
    sns.move_legend(plt.gca(), loc='upper left', bbox_to_anchor=(1,1))

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import scgenome
from scgenome.tools.ranges import dataframe_to_pyranges


def load_snvs(result_adata, vartrix_path):
    vartrix = pd.read_csv(vartrix_path, dtype={'chromosome':str})
    snv_pr = vartrix[['chromosome', 'position']].drop_duplicates().rename(columns={'chromosome':'chr', 'position':'start'})
    snv_pr['end'] = snv_pr['start'] + 1
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        snv_pr = dataframe_to_pyranges(snv_pr)
        region_pr = dataframe_to_pyranges(result_adata.var[['chr', 'start', 'end']])
    amplified_snvs = pyranges_to_dataframe(snv_pr.intersect(region_pr)).rename(columns={'chr':'chromosome', 'start':'position'})[['chromosome', 'position']]
    
    x = vartrix.merge(amplified_snvs, on=['chromosome', 'position'], how='inner')
    print(len(vartrix), len(amplified_snvs), len(x))
    
    alts = x.pivot(index='cell_id', columns=['chromosome', 'position', 'ref', 'alt'], values = 'alt_count')
    refs = x.pivot(index='cell_id', columns=['chromosome', 'position', 'ref', 'alt'], values = 'ref_count')
    
    obs = pd.DataFrame(index=alts.index)
    var = pd.DataFrame(index=alts.columns)
    var = var.reset_index()
    var['snv_id'] = var['chromosome'].astype(str) + ':' + var['position'].astype(str) + ':' + var['ref'] + '>' + var['alt']
    var = var.set_index('snv_id')
    snv_adata = ad.AnnData(obs=obs, var=var, X=alts.values)
    snv_adata.layers['alt'] = alts.loc[snv_adata.obs.index, [tuple(r[['chromosome', 'position', 'ref', 'alt']]) for _, r in snv_adata.var.iterrows()]].fillna(0).values
    snv_adata.layers['ref'] = refs.loc[snv_adata.obs.index, [tuple(r[['chromosome', 'position', 'ref', 'alt']]) for _, r in snv_adata.var.iterrows()]].fillna(0).values
    snv_adata.layers['total'] = snv_adata.layers['alt'] + snv_adata.layers['ref']
    snv_adata.layers['VAF'] = snv_adata.layers['alt'] / snv_adata.layers['total']
    orig_snv_adata = snv_adata.copy()
    
    snv_adata = snv_adata[result_adata.obs.index].copy()
    snv_adata.var['total_counts'] = np.nansum(snv_adata.layers['total'], axis=0)
    snv_adata.var['total_alt_counts'] = np.nansum(snv_adata.layers['alt'], axis=0)
    snv_adata.var['n_cells_without'] = np.sum(np.nan_to_num(snv_adata.layers['VAF'], nan=-1) == 0, axis=0)
    snv_adata.var['n_cells_with'] = np.sum(np.nan_to_num(snv_adata.layers['VAF'], nan=-1) > 0.5, axis=0)
    snv_adata.var['n_cells_missing'] = np.sum(np.isnan(snv_adata.layers['VAF']), axis=0)
    
    contains_snv = np.zeros(snv_adata.shape)
    contains_snv[np.where(np.nan_to_num(snv_adata.layers['VAF'], nan=-1) == 0)] = -1
    contains_snv[np.where(np.nan_to_num(snv_adata.layers['VAF'], nan=-1) > 0.5)] = 1
    snv_adata.layers['contains_snv'] = contains_snv
    return snv_adata

def plot_amplicons_pretty(result_adata, realval_amplicons=True, ax=None, legend=True, include_shading=True,
marker='.', s=None):
    if ax is None:
        ax = plt.gca()
        
    amplicons = result_adata.varm['amplicons'].T
    if realval_amplicons:
        xs = np.arange(amplicons.shape[1])
        for i in range(amplicons.shape[0]):
            nz = np.where(amplicons[i] > 0)[0]
            ax.plot(xs[nz], i + amplicons[i][nz] / (np.max(amplicons)+1), c=plt.get_cmap('tab10')(i/10), label=f'amplicon {i}', linewidth=5)
        ax.set_xlabel('bin')
        ax.set_ylabel('Number of copies')
    else:
        for i in range(amplicons.shape[0]):
            nz = np.where(amplicons[i] > 0)[0]
            ax.scatter(nz, [i] * len(nz), color=plt.get_cmap('tab10')(i/10), label=f'amplicon {i}', marker=marker, s=s)
        ax.set_xlabel('10kb bin')
        ax.set_ylabel('Amplicon index')

    _, ymax = ax.get_ylim()
    boundaries = []
    for i, (amplicon_id, avar) in enumerate(result_adata.var.groupby('amplicon_id', observed=True)):
        start = avar.iindex.min() - 0.5
        end = avar.iindex.max() + 0.5
        if include_shading:
            ax.add_patch(Rectangle((start, 0), 
                            width=(end-start), 
                            height=ymax, zorder=-1, label=amplicon_id, facecolor = plt.get_cmap('tab20')(i), alpha=0.5))
            for discontinuity in np.where(np.diff(avar.iindex) != 1)[0]:
                ax.axvline(x=discontinuity + start + 1, color='grey', linestyle='--')
        boundaries.append(end)
        
    if legend:
        ax.legend()
        sns.move_legend(ax, loc='upper left', bbox_to_anchor=(1,1))
    return boundaries


def plot_amplicon_figure(result_adata, axs, cn_vmax=None, loading_vmax = None, amplicon_marker='|', amplicon_markersize=100,
                        binary_amplicon_plot=False, cn_cmap='plasma', log_scale=False):
    if cn_vmax is None:
        cn_vmax = np.nanmax(result_adata.layers['data'])
   
    if loading_vmax is None:
        loading_vmax = np.nanmax(result_adata.obsm['loadings']) 
    
    n_amplicons = result_adata.obsm['loadings'].shape[1]
    result_adata.var['chr_code'] = result_adata.var['chr'].cat.codes
    ordered_amplicons = result_adata.var.groupby('amplicon_id', observed=True)[['chr_code', 'start']].min().sort_values(by=['chr_code', 'start'])
    n_regions = len(ordered_amplicons)
    binary_amplicons = np.all(result_adata.varm['amplicons'] <= 1)

    if log_scale:
        result_adata.layers['data'] = np.log10(result_adata.layers['data'] + 1)
        print('old vmax:', cn_vmax)
        cn_vmax=np.log10(cn_vmax + 1)
        print('new vmax:', cn_vmax)

    for i, aid in enumerate(ordered_amplicons.index):
        svar = result_adata.var[result_adata.var['amplicon_id'] == aid]
        if binary_amplicons or binary_amplicon_plot:
            plot_amplicons(result_adata[result_adata.obs.sort_values('archetype').index, svar.index], realval_amplicons=False,
                                                ax=axs[1][i], include_shading=False, marker=amplicon_marker, s=amplicon_markersize)
            axs[1][i].set_ylim(-0.5, n_amplicons - 0.5)
        else:
            plot_amplicons_pretty(result_adata[result_adata.obs.sort_values('archetype').index, svar.index], realval_amplicons=True,
                                                ax=axs[1][i], include_shading=False, marker=amplicon_marker, s=amplicon_markersize)
            max_mult = np.max(result_adata.varm['amplicons'])
            axs[1][i].set_ylim(0, n_amplicons)
        axs[1][i].set_xlim(0, len(svar))
        axs[1][i].get_legend().remove()
        axs[1][i].set_xticks([], [])
        axs[1][i].set_yticks([], [])
        axs[1][i].set_xlabel('')
        axs[1][i].set_ylabel('')
        for spine in ['top', 'right', 'left', 'bottom']:
            axs[1][i].spines[spine].set_visible(False)
            
        sns.heatmap(result_adata[result_adata.obs.sort_values('archetype').index, svar.index].layers['data'], ax=axs[0][i], 
                    vmin=0, vmax=cn_vmax, cmap=cn_cmap, rasterized=True, cbar_ax=axs[1][n_regions],
                    cbar_kws={'orientation': 'horizontal', 'pad':1, 'shrink':0.5, 'label':'log10(Copy number + 1)' if log_scale else 'Copy number'})
        if log_scale:
            print(axs[1][n_regions].get_xticks())
            print(axs[1][n_regions].get_xticklabels())
            print(axs[1][n_regions].set_title('AAAA'))
            #axs[1][n_regions].set_xticklabels([f'{int(10**x)}' if x > 1 else '0' for x in axs[1][n_regions].get_xticks()])
        axs[0][i].set_xticks([], [])
        axs[0][i].set_yticks([], [])
        axs[0][i].set_xlabel('')
        axs[0][i].set_ylabel('')

        for spine in ['top', 'right', 'left']:
            axs[2][i].spines[spine].set_visible(False)

        region_start = svar['start'].min() - 1
        region_end = svar['end'].max()
        axs[2][i].set_yticks([], [])
        axs[2][i].set_xlim(region_start, region_end)
        
        # minor tick every 10 bins (100kb)
        minor_checkpoints = []
        minor_checkpoints += [x for x in np.arange(region_start - (region_start % 1e5), region_end, 1e5) if x > region_start]
        axs[2][i].set_xticks(minor_checkpoints, [], minor=True)

        # major tick at the start and every 500kb
        major_checkpoints = [region_start]
        major_checkpoints += [x for x in np.arange(region_start - (region_start % 5e5), region_end, 5e5) if x > region_start]
        axs[2][i].set_xticks(major_checkpoints, [(f'{int(x/1e6)}M' if x % 1e6 == 0 else f'{x/1e6}M') for x in major_checkpoints])

        chrom = svar["chr"].iloc[0]
        axs[2][i].set_xlabel(f'chr{chrom}')
        
        for gene_name, (gene_start, gene_end) in gene_bounds[chrom]:
            axs[2][i].add_patch(Rectangle(xy=(gene_start, 0), width=(gene_end - gene_start), height=0.8))
            axs[2][i].annotate(gene_name, ((gene_start + gene_end)/2, 1), ha='center', va='bottom', style='italic')
        axs[2][i].set_ylim(0, 1)
    axs[1][n_regions].set_ylim(0, 3)
    #axs[1][len(ordered_amplicons.index)].annotate('Copy number', (cn_vmax/2, 1), ha='center', va='bottom')
        

    if binary_amplicons or binary_amplicon_plot:
        axs[1][0].set_yticks(np.arange(1,n_amplicons + 1) - 1, np.arange(1,n_amplicons + 1))
    else:    
        axs[1][0].set_yticks(np.arange(1,n_amplicons + 1) - 0.5, np.arange(1,n_amplicons + 1))

    # first, plot loading heatmap in grayscale for legend
    # plotting this first so that the label for CN colorbar plotted above doesn't get clipped
    loadings = result_adata[result_adata.obs.sort_values('archetype').index].obsm['loadings'].copy()
    sns.heatmap(data=loadings,
                cmap='Greys', vmin=0, vmax=loading_vmax, ax=axs[0][n_regions], rasterized=True,
                cbar_ax=axs[1][n_regions + 1],  cbar_kws={'orientation': 'horizontal', 'label':'Loading'})
    axs[1][n_regions + 1].set_ylim(0, 3)
    axs[1][n_regions + 1].set_zorder(-1)


    for col in np.arange(loadings.shape[1]):
        cmap = mcolors.LinearSegmentedColormap.from_list('', [[0, 'white'], [1, plt.get_cmap('tab10')(col/10)]])
        mask = loadings.copy()
        for j in np.arange(loadings.shape[1]):
            mask[:, j] = col != j
    
        sns.heatmap(data=loadings,
                    mask=mask,
                    cmap=cmap, vmin=0, vmax=loading_vmax, cbar=False, ax=axs[0][n_regions], rasterized=True)


    axs[0][n_regions].set_xticks([], [])
    axs[0][n_regions].set_yticks([], [])

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.2, wspace=0.15)

def get_axes(n_amplicons, width_ratios, dpi=300, figsize=(12,6), height_ratios=[0.9, 0.075, 0.075, 0.02]):
    assert len(height_ratios) == 4, "height_ratios must be a list of 4 elements"
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ncols = len(width_ratios)
    gs = fig.add_gridspec(ncols=ncols, nrows=4, width_ratios = width_ratios,
                           height_ratios = height_ratios)
    
    # assemble array of axes like before to preserve interface
    axs = []
    # first row (main data) is normal
    axs.append([fig.add_subplot(gs[0, j]) for j in range(ncols)])
    
    # second row (amplicon plots and legends)
    # merge rows 1 and 2 for the amplicon axes
    row = [fig.add_subplot(gs[1:3, j]) for j in range(n_amplicons)]
    # split them for the colorbar axes
    row.append(fig.add_subplot(gs[1, n_amplicons]))
    row.append(fig.add_subplot(gs[2, n_amplicons]))

    for i in range(n_amplicons + 1, ncols):
        row.append(fig.add_subplot(gs[1, i]))
    axs.append(row)
    
    # last row is normal, only need an axis for each region
    axs.append([fig.add_subplot(gs[3, j]) for j in range(n_amplicons)])
    return axs
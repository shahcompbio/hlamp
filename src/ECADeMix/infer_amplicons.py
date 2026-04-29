import os
import scgenome
import numpy as np
import pandas as pd
import anndata as ad
import click

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def is_none(x):
    return x is None or x == '' or x == 'nan' or x == 'None'

def load_archetypes(result_adata, patient_id):
    # load manually curated cell orderings for plotting
    loadings = result_adata.obsm['loadings'].copy()
    if patient_id == 'PM_0510':
        # first, group by EGFR amplicon
        archetypes = np.ones(result_adata.shape[0]) * 3000
        archetypes[np.where(loadings[:, 3] > 1)[0]] = 2000
        archetypes[np.where(loadings[:, 1] > 1)[0]] = 1000
        archetypes[np.where(loadings[:, 2] > 7.5)[0]] = 3000
        archetypes[np.where(np.logical_and(loadings[:, 1] > 10, loadings[:, 2] > 5))[0]] = 4000

        # within major EGFR amplicon group, sort by CDK4 +/- and TERT +/-
        archetypes[np.where(loadings[:, 4] < 2)[0]] -= 10
        archetypes[np.where(loadings[:, 0] < 2)[0]] -= 20

        # then, order by TERT and CDK4 amplicons
        archetypes += loadings[:, 4] / np.max(loadings[:, 4])
        archetypes += 0.1 * loadings[:, 0] / np.max(loadings[:, 0])

        result_adata.obs['archetype'] = archetypes

    elif patient_id == 'NCI-H69':
        presence_threshold = 3

        archetypes = np.ones(result_adata.shape[0]) * 500

        single_present = np.where(np.sum(result_adata.obsm['loadings'] > presence_threshold, axis=1) == 1)[0]
        archetypes[single_present] = 100

        double_present = np.where(np.sum(result_adata.obsm['loadings'] > presence_threshold, axis=1) == 2)[0]
        archetypes[double_present] = 200


        archetypes[result_adata.obsm['loadings'][:, 0] > presence_threshold] += 1
        archetypes[result_adata.obsm['loadings'][:, 1] > presence_threshold] += 5
        archetypes[result_adata.obsm['loadings'][:, 2] > presence_threshold] += 10
        archetypes[result_adata.obsm['loadings'][:, 3] > presence_threshold] += 20

        archetypes += np.sum(result_adata.obsm['loadings'], axis=1) / np.max(np.sum(result_adata.obsm['loadings'], axis=1))

        result_adata.obs['archetype'] = archetypes

    elif patient_id == 'Lx298':
        archetypes = np.ones(result_adata.shape[0]) * 300
        archetypes[result_adata.obsm['loadings'][:, 1] < 5] = 600
        archetypes += result_adata.obsm['loadings'][:, 1]
        archetypes += result_adata.obsm['loadings'][:, 0]

        result_adata.obs['archetype'] = archetypes


    elif patient_id == 'NCI-H524':
        cold = np.sum(result_adata.obsm['loadings'] > 5, axis=1) == 0
        onehot = np.sum(result_adata.obsm['loadings'] > 5, axis=1) == 1
        twohot = np.sum(result_adata.obsm['loadings'] > 5, axis=1) == 2

        archetypes = np.ones(result_adata.shape[0]) * 600
        archetypes[cold] = 0

        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 0] > 5)] = 100
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 1] > 5)] = 200
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 2] > 5)] = 300
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 3] > 5)] = 400

        archetypes[np.logical_and(twohot, result_adata.obsm['loadings'][:, 0] > 5)] = 700
        archetypes[np.logical_and(twohot, result_adata.obsm['loadings'][:, 1] > 5)] = 800
        archetypes[np.logical_and(twohot, result_adata.obsm['loadings'][:, 2] > 5)] = 900
        archetypes[np.logical_and(twohot, result_adata.obsm['loadings'][:, 3] > 5)] = 1000

        archetypes += result_adata.obsm['loadings'][:, 0]/np.max(result_adata.obsm['loadings'][:, 0])
        archetypes += result_adata.obsm['loadings'][:, 1]/np.max(result_adata.obsm['loadings'][:, 1])
        archetypes += result_adata.obsm['loadings'][:, 2]/np.max(result_adata.obsm['loadings'][:, 2])
        archetypes += result_adata.obsm['loadings'][:, 3]/np.max(result_adata.obsm['loadings'][:, 3])

        result_adata.obs['archetype'] = archetypes

    elif patient_id == '20200721':
        archetypes = np.ones(result_adata.shape[0]) * 500
        archetypes[(result_adata.obsm['loadings'][:, 1] > 0)] = 100
        archetypes[(result_adata.obsm['loadings'][:, 2] > 0)] = 200
        archetypes[(result_adata.obsm['loadings'][:, 3] > 0)] = 300
        archetypes += result_adata.obsm['loadings'][:, 0]/result_adata.obsm['loadings'][:, 0].max()
        result_adata.obs['archetype'] = archetypes

    elif patient_id == 'P-0009535':
        archetypes = np.ones(result_adata.shape[0]) * 300
        archetypes[result_adata.obsm['loadings'][:, 1] < 1] = 200
        archetypes += result_adata.obsm['loadings'][:, 0]/100
        archetypes += result_adata.obsm['loadings'][:, 1]/5
        result_adata.obs['archetype'] = archetypes

    elif patient_id == 'GBM39-DM':
        archetypes = (result_adata.obsm['loadings'][:, 0] > 0) * 1000 + result_adata.obsm['loadings'][:, 1]
        result_adata.obs['archetype'] = archetypes

    elif patient_id == '2765_2':
        archetypes = np.ones(result_adata.obsm['loadings'].shape[0]) * 500
        archetypes[np.where(result_adata.obsm['loadings'][:, 0] > 5)[0]] = 100
        archetypes[np.where(result_adata.obsm['loadings'][:, 2] > 5)[0]] = 200
        archetypes += result_adata.obsm['loadings'][:, 0]/np.max(result_adata.obsm['loadings'][:, 0])
        archetypes += result_adata.obsm['loadings'][:, 1]/np.max(result_adata.obsm['loadings'][:, 1])
        archetypes += result_adata.obsm['loadings'][:, 2]/np.max(result_adata.obsm['loadings'][:, 2])

        result_adata.obs['archetype'] = archetypes


    elif patient_id == 'COLO320DM':
        cold = np.sum(result_adata.obsm['loadings'] > 10, axis=1) == 0
        onehot = np.sum(result_adata.obsm['loadings'] > 10, axis=1) == 1
        twohot = np.sum(result_adata.obsm['loadings'] > 10, axis=1) == 2

        archetypes = np.ones(result_adata.shape[0]) * 2000
        archetypes[cold] = 0

        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 0] > 10)] = 100
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 1] > 10)] = 200
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 2] > 10)] = 300
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 3] > 10)] = 400
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 4] > 10)] = 500
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 5] > 10)] = 600
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 6] > 10)] = 700
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 7] > 10)] = 800
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 8] > 10)] = 900
        archetypes[np.logical_and(onehot, result_adata.obsm['loadings'][:, 9] > 10)] = 1000

        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 0] > 1, result_adata.obsm['loadings'][:, 1] > 1))] = 1010
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 0] > 1, result_adata.obsm['loadings'][:, 2] > 1))] = 1020
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 0] > 1, result_adata.obsm['loadings'][:, 3] > 1))] = 1030
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 0] > 1, result_adata.obsm['loadings'][:, 4] > 1))] = 1040
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 0] > 1, result_adata.obsm['loadings'][:, 5] > 1))] = 1050
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 0] > 1, result_adata.obsm['loadings'][:, 6] > 1))] = 1060
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 0] > 1, result_adata.obsm['loadings'][:, 7] > 1))] = 1070
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 0] > 1, result_adata.obsm['loadings'][:, 8] > 1))] = 1080

        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 1] > 1, result_adata.obsm['loadings'][:, 2] > 1))] = 1120
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 1] > 1, result_adata.obsm['loadings'][:, 3] > 1))] = 1130
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 1] > 1, result_adata.obsm['loadings'][:, 4] > 1))] = 1140
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 1] > 1, result_adata.obsm['loadings'][:, 5] > 1))] = 1150
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 1] > 1, result_adata.obsm['loadings'][:, 6] > 1))] = 1160
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 1] > 1, result_adata.obsm['loadings'][:, 7] > 1))] = 1170
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 1] > 1, result_adata.obsm['loadings'][:, 8] > 1))] = 1180

        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 2] > 1, result_adata.obsm['loadings'][:, 3] > 1))] = 1230
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 2] > 1, result_adata.obsm['loadings'][:, 4] > 1))] = 1240
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 2] > 1, result_adata.obsm['loadings'][:, 5] > 1))] = 1250
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 2] > 1, result_adata.obsm['loadings'][:, 6] > 1))] = 1260
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 2] > 1, result_adata.obsm['loadings'][:, 7] > 1))] = 1270
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 2] > 1, result_adata.obsm['loadings'][:, 8] > 1))] = 1280

        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 3] > 1, result_adata.obsm['loadings'][:, 4] > 1))] = 1340
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 3] > 1, result_adata.obsm['loadings'][:, 5] > 1))] = 1350
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 3] > 1, result_adata.obsm['loadings'][:, 6] > 1))] = 1360
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 3] > 1, result_adata.obsm['loadings'][:, 7] > 1))] = 1370
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 3] > 1, result_adata.obsm['loadings'][:, 8] > 1))] = 1380

        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 4] > 1, result_adata.obsm['loadings'][:, 5] > 1))] = 1450
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 4] > 1, result_adata.obsm['loadings'][:, 6] > 1))] = 1460
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 4] > 1, result_adata.obsm['loadings'][:, 7] > 1))] = 1470
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 4] > 1, result_adata.obsm['loadings'][:, 8] > 1))] = 1480

        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 5] > 1, result_adata.obsm['loadings'][:, 6] > 1))] = 1560
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 5] > 1, result_adata.obsm['loadings'][:, 7] > 1))] = 1570
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 5] > 1, result_adata.obsm['loadings'][:, 8] > 1))] = 1580

        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 6] > 1, result_adata.obsm['loadings'][:, 6] > 1))] = 1660
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 6] > 1, result_adata.obsm['loadings'][:, 7] > 1))] = 1670
        archetypes[np.logical_and(twohot, np.logical_and(result_adata.obsm['loadings'][:, 6] > 1, result_adata.obsm['loadings'][:, 8] > 1))] = 1680

        archetypes += result_adata.obsm['loadings'][:, 0]/np.max(result_adata.obsm['loadings'][:, 0])
        archetypes += (10**-1) * result_adata.obsm['loadings'][:, 1]/np.max(result_adata.obsm['loadings'][:, 1])
        archetypes += (10**-2) * result_adata.obsm['loadings'][:, 2]/np.max(result_adata.obsm['loadings'][:, 2])
        archetypes += (10**-3) * result_adata.obsm['loadings'][:, 3]/np.max(result_adata.obsm['loadings'][:, 3])
        archetypes += (10**-4) * result_adata.obsm['loadings'][:, 4]/np.max(result_adata.obsm['loadings'][:, 4])
        archetypes += (10**-5) * result_adata.obsm['loadings'][:, 5]/np.max(result_adata.obsm['loadings'][:, 5])
        archetypes += (10**-6) * result_adata.obsm['loadings'][:, 6]/np.max(result_adata.obsm['loadings'][:, 6])
        archetypes += (10**-7) * result_adata.obsm['loadings'][:, 7]/np.max(result_adata.obsm['loadings'][:, 7])
        archetypes += (10**-8) * result_adata.obsm['loadings'][:, 8]/np.max(result_adata.obsm['loadings'][:, 8])
        archetypes += (10**-9) * result_adata.obsm['loadings'][:, 9]/np.max(result_adata.obsm['loadings'][:, 9])
        archetypes += (10**-10) * result_adata.obsm['loadings'][:, 10]/np.max(result_adata.obsm['loadings'][:, 10])

        result_adata.obs['archetype'] = archetypes

    elif patient_id in ['Lx516', 'Lx33', 'P-0093087'] or result_adata.obsm['loadings'].shape[1] == 1:
        # patients with simple 1-amplicon solutions
        archetypes = result_adata.obsm['loadings'][:, 0]
        result_adata.obs['archetype'] = archetypes

    else:
        print(f'no archetypes for patient {patient_id}')
        result_adata.obs['archetype'] = 0

    return result_adata

@click.command()
@click.option("--patient_id", required=True)
@click.option("--ilp_input_adata", required=True)
@click.option("--hmmcopy_reads", required=True, multiple=True)
@click.option("--hmmcopy_metrics", required=True, multiple=True)
@click.option("--patient_id", required=True)
@click.option("--min_quality", required=True, type=float)
@click.option("--max_reads_nonamplicon_cell", required=True, type=float)
@click.option("--normalize_data", is_flag=True, default=False)
@click.option("--round_cn_normalize", is_flag=True, default=False)
@click.option("--noamplicon_profile_file", required=False, default=None)
@click.option("--manual_init_file", required=False, default=None)
@click.option("--n_amplicons", required=True, type=int)
@click.option("--max_amplicon_cn", required=True, type=int)
@click.option("--max_loading", required=True, type=int)
@click.option("--min_loading", required=True, type=int)
@click.option("--lambda_x", required=True, type=float)
@click.option("--lambda_x_l1", required=True, type=float)
@click.option("--lambda_y", required=True, type=float)
@click.option("--cell_weights", is_flag=True, default=False)
@click.option("--refit_loadings_lambda_y", required=False, type=float, default=np.nan)
@click.option("--plot_data_vmax", required=False, default=500, type=float)
@click.option("--plot_loading_vmax", required=False, default=100, type=float)
@click.option("--seed", required=True, type=int)
@click.option("--gurobi_license_file", required=True)
@click.option("--outdir", required=True)
def construct_clones(
    patient_id,
    ilp_input_adata,
    hmmcopy_reads,
    hmmcopy_metrics,
    min_quality ,
    max_reads_nonamplicon_cell,
    normalize_data,
    round_cn_normalize,
    noamplicon_profile_file,
    manual_init_file,
    n_amplicons,
    max_amplicon_cn,
    max_loading,
    min_loading,
    lambda_x,
    lambda_x_l1,
    lambda_y,
    cell_weights,
    refit_loadings_lambda_y,
    plot_data_vmax,
    plot_loading_vmax,
    seed,
    gurobi_license_file,
    outdir,     
):
    ### sanitize inputs
    if is_none(noamplicon_profile_file):
        noamplicon_profile_file = None
    if is_none(manual_init_file):
        manual_init_file = None
    if is_none(cell_weights):
        cell_weights = False
    if is_none(round_cn_normalize):
        round_cn_normalize = False
    if is_none(normalize_data):
        normalize_data = False

    ###
    print('normalize data: ', normalize_data)
    print('round CN normalize: ', round_cn_normalize)
    
    # set Gurobi license before import GRB packages (imported by ecademix)
    os.environ['GRB_LICENSE_FILE'] = gurobi_license_file # must precede ecademix import
    import ecademix.solve
    import ecademix.postprocess

    # load patient anndata
    adata = ad.read_h5ad(ilp_input_adata)

    # restrict to amplicons, sort and annotate
    amplicon_adata = adata[:, adata.var.in_amplicon].copy()
    amplicon_adata = amplicon_adata[:, amplicon_adata.var.sort_values(by=['amplicon_id', 'start']).index]
    amplicon_adata.var = amplicon_adata.var.reset_index(drop=True)
    amplicon_adata.var['iindex'] = np.arange(len(amplicon_adata.var))
    amplicon_adata.var['total_counts'] = np.nansum(amplicon_adata.layers['reads'], axis=0)
    amplicon_adata.var['total_cpm'] = np.nansum(amplicon_adata.layers['cpm'], axis=0)



    # load 500kb hmmcopy and restrict to high-quality cells
    metrics = pd.concat([pd.read_csv(f) for f in hmmcopy_metrics])
    reads = pd.concat([pd.read_csv(f) for f in hmmcopy_reads])
    hmmcopy = scgenome.pp.convert_dlp_hmmcopy(metrics, reads)
    hmmcopy.obs['ploidy'] = np.nanmean(hmmcopy.layers['state'], axis=1)
    hq = hmmcopy[(hmmcopy.obs.quality > min_quality) & (hmmcopy.obs.index.isin(amplicon_adata.obs.index))].copy()
    rpc = []
    for cell in hq.obs.index:
        valid_idx = np.where(hq[cell].layers['state'][0] > 0)[0]
        rpc.append(np.mean(hq[cell].X[0][valid_idx] / hq[cell].layers['state'][0][valid_idx]))
    hq.obs['mean_reads_per_copy'] = rpc
    hq.obs = hq.obs.merge(amplicon_adata.obs, left_index=True, right_index=True)

    # scale 10kb reads to copies using proportionality from 500kb profile
    data2 = ((amplicon_adata[hq.obs.index].layers['reads'].T * (5e5/1e4)) / hq.obs['mean_reads_per_copy'].values).T
    hq.obs['row_sum'] = row_sum = np.sum(data2, axis=1)

    print(f'number of NaNs in data2: {np.sum(np.isnan(data2))}')


    # normalize using amplicon-free cells if present
    if normalize_data:
        normal_profile = np.nanmean(data2[row_sum <= max_reads_nonamplicon_cell], axis=0)
        orig_profile = normal_profile.copy()
        if round_cn_normalize:
            normal_profile /= np.maximum(1, np.round(normal_profile))
        if noamplicon_profile_file:
            nonamplicon_profile = np.loadtxt(noamplicon_profile_file, delimiter=',')
            normal_profile /= nonamplicon_profile
        normal_profile /= np.mean(normal_profile)
        data3 = data2 / normal_profile

        # plot data before and after normaliation?
        plt.figure(figsize=(8,2), dpi=300)
        plt.plot(orig_profile, label='unnormalized divisor')

        plt.plot(normal_profile, label='normalized divisor')
        _, ymax = plt.ylim()
        for i, (amplicon_id, avar) in enumerate(amplicon_adata.var.groupby('amplicon_id', observed=True)):
            start = avar.iindex.min() - 0.5
            end = avar.iindex.max() + 0.5
            plt.gca().add_patch(Rectangle((start, 0), 
                            width=(end-start), 
                            height=ymax, zorder=-1, label=amplicon_id, facecolor = plt.get_cmap('tab20')(i), alpha=0.5))
            for discontinuity in np.where(np.diff(avar.iindex) != 1)[0]:
                plt.gca().axvline(x=discontinuity + start + 1, color='grey', linestyle='--')
        plt.legend()
        plt.savefig(f'{outdir}/{patient_id}_normalization.png', dpi=300)
        plt.close()
    else:
        data3 = data2.copy()

    # run ILP
    if manual_init_file:
        amplicon_start = np.loadtxt(manual_init_file, delimiter=',')
        loading_start = np.ones((data3.shape[0], n_amplicons))
    else:
        W, H = ecademix.solve.run_NMF(data3, n_amplicons, seed=seed)
        loading_start, amplicon_start = ecademix.solve.construct_solution(W, H, max_amplicon_cn)

    print(f'n_amplicons: {n_amplicons}')
    print(f'max_amplicon_cn: {max_amplicon_cn}')
    print(f'min_loading: {min_loading}')
    print(f'max_loading: {max_loading}')
    print(f'lambda_x: {lambda_x}')
    print(f'lambda_x_l1: {lambda_x_l1}')
    print(f'lambda_y: {lambda_y}')
    print(f'cell_weights: {cell_weights}')
    print(f'type(cell_weights): {type(cell_weights)}')
    print(f'refit_loadings_lambda_y: {refit_loadings_lambda_y}')
    print(f'seed: {seed}')

    result = ecademix.solve.optimize_model(data=data3, 
                                        loading_start=loading_start, 
                                        amplicon_start=amplicon_start, 
                                        n_amplicons=n_amplicons, 
                                        min_loading=min_loading, 
                                        max_loading=max_loading, 
                                        max_amplicon_cn=max_amplicon_cn,
                                        lambda_x=lambda_x, 
                                        lambda_x_l1=lambda_x_l1,
                                        lambda_y=lambda_y, 
                                        verbose=True,
                                        cell_weights=cell_weights)

    result_adata = ecademix.postprocess.compose_results(hq, amplicon_adata, data3, result)

    result_adata = load_archetypes(result_adata, patient_id)

    result_adata.write_h5ad(f'{outdir}/{patient_id}_result_adata.h5ad')

    # QC plots
    ecademix.postprocess.plot_residuals(result_adata, realval_amplicons=True, data_vmax=plot_data_vmax, loading_vmax=plot_loading_vmax)

if __name__ == "__main__":
    construct_clones()

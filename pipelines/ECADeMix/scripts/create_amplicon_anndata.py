import numpy as np
import pandas as pd
import anndata as ad
import click

@click.command()
@click.option("--patient_id", required=True)
@click.option("--regions_table", required=True)
@click.option("--counts_anndata", required=True)
@click.option("--amplicon_anndata", required=True)   
def create_amplicon_anndata(
    patient_id,
    regions_table,
    counts_anndata,
    amplicon_anndata,    
):

    regions = pd.read_table(regions_table)
    adata = ad.read_h5ad(counts_anndata)
    
    adata.var['chr'] = [a.split(':')[0] for a in adata.var.index]
    adata.var['start'] = [a.split(':')[1].split('-')[0] for a in adata.var.index]
    adata.var['end'] = [a.split(':')[1].split('-')[1] for a in adata.var.index]
    adata.var.start = adata.var.start.astype(int)
    adata.var.end = adata.var.end.astype(int)
    adata.layers['reads'] = adata.X
        
    adata.var['in_amplicon'] = False
    adata.var['amplicon_id'] = ''
    my_regions = regions[regions['individual'] == patient_id]
    for _, r in my_regions.iterrows():
        # include flanking regions of 1mb on each side for 20200721
        region_start = r.start - (r.start % 10000) - (1e6 if patient_id == '20200721' else 0)
        region_end = r.end + (10000 - r.end % 10000) + (1e6 if patient_id == '20200721' else 0)
        region_key = (adata.var.chr == str(r.chromosome)) & (adata.var.start >= region_start) & (adata.var.end <= region_end)
        adata.var.loc[region_key, 'in_amplicon'] = True
        adata.var.loc[region_key, 'amplicon_id'] =  f'{r.chromosome}:{r.start}-{r.end}'
        print(f'{patient_id} chr{r.chromosome} actual={r.end - r.start}, rounded={region_end - region_start}, {region_key.sum()} 10kb bins')
    
    adata.obs['total_counts'] = np.sum(adata.X, axis = 1).astype(int)
    adata.layers['cpm'] = (1e6 * (adata.X.T / np.maximum(1, adata.obs.total_counts.values)).T).astype(float)
    adata.var['mean_cpm'] = np.mean(adata.layers['cpm'], axis = 0).astype(float)
    adata.var['variance'] = np.var(adata.X, axis = 0)

    
    adata = adata[:, adata.var[adata.var.in_amplicon].index]
    
    if patient_id == '20200721':
        # restrict flanking regions to only those that are amplified in some cells (determined by manual inspection)
        adata = adata[:, np.concatenate([np.arange(99, 217), np.arange(370, 525)])].copy()
    
    adata.obs['amplicon_total_reads'] = np.sum(adata.layers['reads'], axis=1)
    adata.var['total_counts'] = np.nansum(adata.layers['reads'], axis=0)
    adata.var['total_cpm'] = np.nansum(adata.layers['cpm'], axis=0)
    adata = adata[
        (adata.obs.total_counts > 0) & 
        (adata.obs.amplicon_total_reads > 0),
        (adata.var.total_counts > 0) & 
        (adata.var.total_cpm > 0)
    ].copy()

    adata.write_h5ad(amplicon_anndata)

if __name__ == "__main__":
    create_amplicon_anndata()

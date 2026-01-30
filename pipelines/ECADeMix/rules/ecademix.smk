rule get_ilp_input:
    input:
        expand(outdir + '/ecademix/input/{patient_id}_ilp_input.h5ad', patient_id = patients)

rule run_all_ecademix:
    input:
        expand(outdir + '/ecademix/output/{patient_id}/{patient_id}_result_adata.h5ad', patient_id = patients)

rule create_amplicon_anndata:
    input:
        regions_table=config['regions_table'],
        counts_anndata=outdir + '/cell-counts/{patient_id}_cell_counts.h5'
    output:
        ilp_input_adata=outdir + '/ecademix/input/{patient_id}_ilp_input.h5ad'
    resources:
        mem_mb=64 * 1000,
        runtime=59,
        partition='cpushort',
        slurm_partition='cpushort',
    conda: 'spectrum'
    shell:
        '''
        {python_bin} {scripts_dir}/create_amplicon_anndata.py \
            --patient_id {wildcards.patient_id} \
            --regions_table {input.regions_table} \
            --counts_anndata {input.counts_anndata} \
            --amplicon_anndata {output.ilp_input_adata}
        '''


def get_init_file(wildcards):
    if ecademix_params.loc[wildcards.patient_id, 'manual_init_file'] is not None:
        return metadata + f'/{wildcards.patient_id}_amplicons.txt'
    else:
        return None

def get_noamplicon_profile_file(wildcards):
    if ecademix_params.loc[wildcards.patient_id, 'manual_init_file'] is not None:
        return metadata + f'/{wildcards.patient_id}_nonamplicon_divisor.txt'
    else:
        return None

def get_ecademix_input(wildcards):
    hmmcopy_reads = []
    hmmcopy_metrics = []
    for aliquot, adf in hmmcopy[hmmcopy.isabl_patient_id == wildcards.patient_id].groupby('isabl_aliquot_id'):
        adf = adf.sort_values(by='analysis_pk', ascending=False)
        hmmcopy_reads.append(adf.iloc[0].reads)
        hmmcopy_metrics.append(adf.iloc[0].metrics)
    if len(hmmcopy_reads) == 0:
        raise ValueError(f'No hmmcopy reads found for patient {wildcards.patient_id}')
    if len(hmmcopy_metrics) == 0:
        raise ValueError(f'No hmmcopy metrics found for patient {wildcards.patient_id}')
    assert len(hmmcopy_reads) == len(hmmcopy_metrics)

    return {
        'ilp_input_adata': outdir + '/ecademix/input/{patient_id}_ilp_input.h5ad',
        'hmmcopy_reads': hmmcopy_reads,
        'hmmcopy_metrics': hmmcopy_metrics,
    }

rule run_ecademix:
    input:
        unpack(get_ecademix_input)
    params:
        # explicitly formed arguments
        hmmcopy_reads_arg=lambda wildcards, input: ' --hmmcopy_reads '.join(input.hmmcopy_reads),
        hmmcopy_metrics_arg=lambda wildcards, input: ' --hmmcopy_metrics '.join(input.hmmcopy_metrics),
        noamplicon_profile_file_arg=lambda wildcards: f'--noamplicon_profile_file {metadata}/init-files/{wildcards.patient_id}_noamplicon_divisor.txt' if ecademix_params.loc[wildcards.patient_id, 'noamplicon_profile_file'] else '',
        manual_init_file_arg=lambda wildcards: f'--manual_init_file {metadata}/init-files/{wildcards.patient_id}_amplicons.txt' if ecademix_params.loc[wildcards.patient_id, 'manual_init_file'] else '',
        # flags
        cell_weights_arg=lambda wildcards: '--cell_weights' if ecademix_params.loc[wildcards.patient_id, 'cell_weights'] else '',
        normalize_data_arg=lambda wildcards: '--normalize_data' if ecademix_params.loc[wildcards.patient_id, 'normalize_data'] else '',
        round_cn_normalize_arg=lambda wildcards: '--round_cn_normalize' if ecademix_params.loc[wildcards.patient_id, 'round_cn_normalize'] else '',
        # standard arguments
        n_amplicons=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'n_amplicons'],
        max_amplicon_cn=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'max_amplicon_cn'],
        max_loading=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'max_loading'],
        min_loading=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'min_loading'],
        min_quality=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'min_quality'],
        max_reads_nonamplicon_cell=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'max_reads_nonamplicon_cell'],
        lambda_x=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'lambda_x'],
        lambda_x_l1=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'lambda_x_l1'],
        lambda_y=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'lambda_y'],
        plot_data_vmax=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'plot_data_vmax'],
        plot_loading_vmax=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'plot_loading_vmax'],
        seed=lambda wildcards: ecademix_params.loc[wildcards.patient_id, 'seed'],
        gurobi_license_file=gurobi_license_file,
        outdir=lambda wildcards: outdir + f'/ecademix/output/{wildcards.patient_id}'
    output:
        outdir + '/ecademix/output/{patient_id}/{patient_id}_result_adata.h5ad'
    resources:
        mem_mb=64 * 1000,
        runtime=3*24*60,
        partition='componc_cpu',
        slurm_partition='componc_cpu',
    conda: 'spectrum'
    shell:
        '''
        mkdir -p {params.outdir}
        {python_bin} {scripts_dir}/infer_amplicons.py \
            --patient_id {wildcards.patient_id} \
            --ilp_input_adata {input.ilp_input_adata} \
            --hmmcopy_reads {params.hmmcopy_reads_arg} \
            --hmmcopy_metrics {params.hmmcopy_metrics_arg} \
            --min_quality {params.min_quality} \
            --max_reads_nonamplicon_cell {params.max_reads_nonamplicon_cell} \
            {params.noamplicon_profile_file_arg} \
            {params.manual_init_file_arg} \
            --n_amplicons {params.n_amplicons} \
            --max_amplicon_cn {params.max_amplicon_cn} \
            --max_loading {params.max_loading} \
            --min_loading {params.min_loading} \
            --lambda_x {params.lambda_x} \
            --lambda_x_l1 {params.lambda_x_l1} \
            --lambda_y {params.lambda_y} \
            {params.cell_weights_arg} \
            {params.normalize_data_arg} \
            {params.round_cn_normalize_arg} \
            --plot_data_vmax {params.plot_data_vmax} \
            --plot_loading_vmax {params.plot_loading_vmax} \
            --seed {params.seed} \
            --gurobi_license_file {params.gurobi_license_file} \
            --outdir {params.outdir}
        '''

rule plot_ecademix_results:
    input:
        outdir + '/ecademix/output/{patient_id}/{patient_id}_result_adata.h5ad'
    output:
        outdir + '/ecademix/output/{patient_id}/{patient_id}_results.svg'
    shell:
        '''
        {python_bin} {scripts_dir}/plot_ecademix_results.py \
            --patient_id {wildcards.patient_id} \
            --adata {input.adata} \
            --output_file {output.output_file}
        '''
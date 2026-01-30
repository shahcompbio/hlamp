
rule get_counts:
    input:
        expand(outdir + '/cell-counts/{patient}_cell_counts.h5', patient=patients)

rule get_control_counts:
    input:
        expand(outdir + '/cell-counts/{patient}_control_cell_counts.h5', patient=patients)


def get_bamfile(wildcards):
    bams = bamfiles.loc[bamfiles.isabl_aliquot_id == wildcards.aliquot, 'result_filepath'].values
    assert len(bams) == 1
    return bams[0]


rule count_reads:
    input:
        get_bamfile
    output:
        temp(outdir + '/cell-counts/aliquot/aliquot_cell_counts_{aliquot}_chrs{chrs}.csv.gz')
    params:
        binsize = binsize,
        min_mapq = min_mapq,
        chromosomes = lambda wildcards: ' --chromosomes '.join(wildcards.chrs.split('-')),
    threads: 16
    conda: 'mondrianutils'
    benchmark:
        os.path.join(benchmark_dir, 'count_reads', 'count_reads.{aliquot}.{chrs}.txt')
    resources:
        mem_mb=32*1000,# lambda wc, input: max(8000, 100*input.size_mb)
        runtime="1h",
        partition='cpushort',
        slurm_partition='cpushort',
    shell:
        """
        hmmcopy_utils readcounter \
            --infile {input} \
            --outdir {output} \
            -w {params.binsize} \
            -m {params.min_mapq} \
            --ncores 16  --tabular \
            --chromosomes {params.chromosomes}
        """

def get_control_bamfile(wildcards):
    bams = bamfiles.loc[bamfiles.isabl_aliquot_id == wildcards.aliquot, 'result_control_filepath'].values
    assert len(bams) == 1
    return bams[0]
rule count_control_reads:
    input:
        get_control_bamfile
    output:
        temp(outdir + '/cell-counts/aliquot/aliquot_control_cell_counts_{aliquot}_chrs{chrs}.csv.gz')
    params:
        binsize = binsize,
        min_mapq = min_mapq,
        chromosomes = lambda wildcards: ' --chromosomes '.join(wildcards.chrs.split('-')),
    threads: 16
    conda: 'mondrianutils'
    benchmark:
        os.path.join(benchmark_dir, 'count_reads', 'count_reads.{aliquot}.{chrs}.txt')
    resources:
        mem_mb=32*1000,# lambda wc, input: max(8000, 100*input.size_mb)
        runtime="1h",
        partition='cpushort'
    shell:
        """
        hmmcopy_utils readcounter \
            --infile {input} \
            --outdir {output} \
            -w {params.binsize} \
            -m {params.min_mapq} \
            --ncores 16  --tabular \
            --chromosomes {params.chromosomes}
        """
def get_control_aliquot_tables(wildcards):
    # use only the most recent aliquot per sample
    my_aliquots = []
    my_df = bamfiles[bamfiles.isabl_patient_id == wildcards.patient]
    for sample_id in my_df.isabl_sample_id.unique():
        sample_df = my_df[my_df.isabl_sample_id == sample_id].sort_values(by='isabl_pk', ascending=False)
        aliquot = sample_df.isabl_aliquot_id.iloc[0]
        if os.path.getsize(sample_df.iloc[0].result_control_filepath) < min_bamfile_size:
            print("Found aliquot that was too small:", aliquot)
            continue
        else:
            my_aliquots.append(aliquot)
    #my_aliquots = bamfiles[bamfiles.isabl_patient_id == wildcards.patient].isabl_aliquot_id.unique()
    assert len(my_aliquots) > 0, f"Found no aliquots for patient: {wildcards.patient}"
    return [outdir + f'/cell-counts/aliquot/aliquot_control_cell_counts_{a}_chrs{"-".join(chrs)}.csv.gz' 
            for a in my_aliquots for chrs in chr_sets]

rule merge_control_counts:
    input:
        tables=get_control_aliquot_tables
    params:
        tables_arg=lambda wildcards, input:' --aliquot_tables '.join(input.tables),
        aliquot_counts_dir = outdir + '/cell-counts/aliquot/'
    resources:
        mem_mb=156 * 1000,
        runtime="1h",
        partition='cpushort'
    output:
        outdir + '/cell-counts/{patient}_control_cell_counts.h5'
    benchmark:
        os.path.join(benchmark_dir, 'merge_counts', 'merge_counts.{patient}.txt')
    shell:
        '''
        python {scripts_dir}/merge_counts.py --aliquot_tables {params.tables_arg} \
            --patient_id {wildcards.patient} --outfile {output} \
            --aliquot_counts_dir {params.aliquot_counts_dir}
        '''


def get_chr_sets(is_mouse):
    if is_mouse:
        chr_sets = [['1'], ['2'], ['X']]
        remaining_chrs = [str(i) for i in range(3, 20)] + ['Y']
        for i in range(9):
            chr_sets.append([remaining_chrs[i], remaining_chrs[-1 * i - 1]])    
    else:
        # split chromosomes into 13 sets to balance load
        chr_sets = [['1'], ['2']]
        remaining_chrs = [str(i) for i in range(3, 8)] + ['X'] +[str(i) for i in range(8, 23)] + ['Y']
        for i in range(11):
            chr_sets.append([remaining_chrs[i], remaining_chrs[-1 * i - 1]])
    return chr_sets

def get_aliquot_tables(wildcards):
    # use only the most recent aliquot per sample
    my_aliquots = []
    my_df = bamfiles[bamfiles.isabl_patient_id == wildcards.patient]
    chr_sets = get_chr_sets(wildcards.patient in mouse_patients)
    for sample_id in my_df.isabl_sample_id.unique():
        sample_df = my_df[my_df.isabl_sample_id == sample_id].sort_values(by='isabl_pk', ascending=False)
        my_aliquots.append(sample_df.isabl_aliquot_id.iloc[0])
    assert len(my_aliquots) > 0, f"Found no aliquots for patient: {wildcards.patient}"
    return [outdir + f'/cell-counts/aliquot/aliquot_cell_counts_{a}_chrs{"-".join(chrs)}.csv.gz' 
            for a in my_aliquots for chrs in chr_sets]


rule merge_counts:
    input:
        tables=get_aliquot_tables
    params:
        tables_arg=lambda wildcards, input:' --aliquot_tables '.join(input.tables),
        aliquot_counts_dir = outdir + '/cell-counts/aliquot/'
    resources:
        mem_mb=156 * 1000,
        runtime="1h",
        partition='cpushort'
    output:
        outdir + '/cell-counts/{patient}_cell_counts.h5'
    benchmark:
        os.path.join(benchmark_dir, 'merge_counts', 'merge_counts.{patient}.txt')
    shell:
        '''
        python {scripts_dir}/merge_counts.py --aliquot_tables {params.tables_arg} \
            --patient_id {wildcards.patient} --outfile {output} \
            --aliquot_counts_dir {params.aliquot_counts_dir}
        '''

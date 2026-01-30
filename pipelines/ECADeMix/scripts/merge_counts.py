import os
import pandas as pd
import anndata as ad
import tqdm
import click


def load_adata(counts_file):
    df = pd.read_csv(counts_file, dtype={"chromosome": str})
    pt = df.pivot(
        index="cell_id", columns=["chromosome", "start", "end"], values="count"
    )

    obs_df = pd.DataFrame(pt.index.rename("cell_id")).set_index("cell_id")
    var_df = (
        pt.columns.to_series()
        .reset_index()[["chromosome", "start", "end"]]
        .rename(columns={"chromosome": "chr"})
    )
    var_df.chr = var_df.chr.astype(str)
    var_df.start = var_df.start.astype(int) + 1
    var_df.end = var_df.end.astype(int)
    var_df.index = (
        var_df.chr + ":" + var_df.start.astype(str) + "-" + var_df.end.astype(str)
    )
    part_adata = ad.AnnData(obs=obs_df, var=var_df, X=pt.values)
    return part_adata


@click.command()
@click.option("--aliquot_tables", required=True, multiple=True)
@click.option("--aliquot_counts_dir", required=True)
@click.option("--patient_id", required=True)
@click.option("--outfile", required=True)
def merge_counts(aliquot_tables, aliquot_counts_dir, patient_id, outfile):
    pair2table = {}
    aliquots = set()
    chr_sets = set()
    for af in aliquot_tables:
        if not os.path.exists(af):
            raise FileNotFoundError(f"Aliquot table {af} does not exist.")
        chr_set = af.split(".")[0].split("_")[-1]
        aliquot = "_".join(af.split(".")[0].split("_")[3:-1])
        
        pair2table[(aliquot, chr_set)] = af
        aliquots.add(aliquot)
        chr_sets.add(chr_set)

    # first combine horizontally within aliquots
    aliquot_adatas = []
    for aliquot in tqdm.tqdm(aliquots):
        for ch in chr_sets:
            if (aliquot, ch) not in pair2table:
                raise ValueError(
                    f"Missing table for aliquot {aliquot} and chr set {ch}."
                )
            elif not os.path.exists(pair2table[(aliquot, ch)]):
                        raise FileNotFoundError(
                            f"Table {pair2table[(aliquot, ch)]} does not exist."
                )

        aliquot_parts = []
        for ch in chr_sets:
            part_adata = load_adata(
                pair2table[(aliquot, ch)]
            )
            aliquot_parts.append(part_adata)
        aliquot_adata = ad.concat(aliquot_parts, axis=1)
        aliquot_adata.obs["aliquot_id"] = aliquot
        aliquot_adatas.append(aliquot_adata)

    # then combine aliquots
    patient_adata = ad.concat(aliquot_adatas, axis=0)
    patient_adata.obs["patient_id"] = patient_id
    patient_adata.var["chr"] = [a.split(":")[0] for a in patient_adata.var.index]
    patient_adata.var["start"] = [
        a.split(":")[1].split("-")[0] for a in patient_adata.var.index
    ]
    patient_adata.var["end"] = [
        a.split(":")[1].split("-")[1] for a in patient_adata.var.index
    ]
    patient_adata.var.start = patient_adata.var.start.astype(int)
    patient_adata.var.end = patient_adata.var.end.astype(int)
    patient_adata.write(outfile)


if __name__ == "__main__":
    merge_counts()

# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

# %%
# %load_ext autoreload
# %autoreload 2

from nxfvars import nxfvars
import scanpy as sc
import numpy as np
import itertools
from tqdm import trange, tqdm
import scipy.sparse
import numpy.testing as npt
from scanpy_helpers.integration import (
    normalize_by_gene_length,
    sanitize_adata,
    validate_adata,
    add_doublet_annotation,
    undo_log_norm,
    remap_gene_symbols,
    drop_duplicated_genes,
    aggregate_duplicate_gene_symbols,
    merge_datasets,
    MANDATORY_COLS,
)
from threadpoolctl import threadpool_limits
from tqdm.contrib.concurrent import process_map
import mygene
from operator import and_
from functools import reduce
import pandas as pd
import anndata
import re
import os

# %%
out_dir = nxfvars.get("artifact_dir", "/local/scratch/sturm/")

# %%
threadpool_limits(int(nxfvars.get("cpus", "8")))

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
dataset_table = pd.read_csv(
    nxfvars.get("samplesheet", "../../tables/samplesheet_scrnaseq_preprocessing.csv")
)
dataset_path_annotated = nxfvars.get(
    "dataset_path_annotated",
    "../../data/20_build_atlas/integrate_datasets/11_seed_annotations/artifacts/",
)
dataset_path = nxfvars.get(
    "dataset_path", "../../data/20_build_atlas/integrate_datasets/02_qc_and_filtering/"
)

# %%
dataset_table

# %%
datasets_annotated = ["Maynard_Bivona_2020_NSCLC", "Lambrechts_2018_LUAD_6653"]

# %%
datasets = {
    dataset_id: sc.read_h5ad(
        f"{dataset_id}.qc.h5ad"
        if dataset_path == "."
        else f"{dataset_path}/{dataset_id}/{dataset_id}.qc.h5ad"
    )
    for dataset_id in tqdm(dataset_table["id"])
}

# %%
# Set cell-types of unannotated datasets to "unknown" for scVI
for dataset_id in datasets:
    datasets[dataset_id].obs["cell_type"] = "unknown"

# %%
for dataset_id in datasets_annotated:
    tmp_adata = sc.read_h5ad(f"{dataset_path_annotated}/{dataset_id}_annotated.h5ad")
    datasets[dataset_id].obs["cell_type"] = tmp_adata.obs["cell_type"]

# %% [markdown]
# ### Dataset-specific filtering and metadata fixes

# %%
datasets["Adams_Kaminski_2020_COPD"].obs["origin"] = "normal"
datasets["Adams_Kaminski_2020_COPD"].obs["sex"] = "nan"
datasets["Adams_Kaminski_2020_COPD"] = datasets["Adams_Kaminski_2020_COPD"][
    datasets["Adams_Kaminski_2020_COPD"].obs["condition"] != "IPF", :
]

# %%
# No modifications necessary for Chen_Zhang

# %%
datasets["Goveia_Carmeliet_2020_NSCLC"] = datasets["Goveia_Carmeliet_2020_NSCLC"][
    datasets["Goveia_Carmeliet_2020_NSCLC"].obs["condition"] != "LLCC"
].copy()
datasets["Goveia_Carmeliet_2020_NSCLC"].obs["sex"] = "nan"

# %%
datasets["Guo_Zhang_2018_NSCLC"] = datasets["Guo_Zhang_2018_NSCLC"][
    datasets["Guo_Zhang_2018_NSCLC"].obs["tissue"] != "blood"
].copy()
# store true raw counts for later export
tmp_raw_counts = datasets["Guo_Zhang_2018_NSCLC"].X
datasets["Guo_Zhang_2018_NSCLC"] = normalize_by_gene_length(
    datasets["Guo_Zhang_2018_NSCLC"]
)
datasets["Guo_Zhang_2018_NSCLC"].layers["raw_counts"] = tmp_raw_counts
datasets["Guo_Zhang_2018_NSCLC"].obs["sex"] = "nan"

# %%
datasets["Habermann_Kropski_2020_pulmonary-fibrosis"].obs["sex"] = [
    {"M": "male", "F": "female", "Unknown": "nan"}[s]
    for s in datasets["Habermann_Kropski_2020_pulmonary-fibrosis"].obs["sex"]
]
datasets["Habermann_Kropski_2020_pulmonary-fibrosis"] = datasets[
    "Habermann_Kropski_2020_pulmonary-fibrosis"
][
    datasets["Habermann_Kropski_2020_pulmonary-fibrosis"].obs["condition"]
    == "healthy_control",
    :,
].copy()

# %%
# No modifications necessary for He_Fan

# %%
# store true raw counts for later export
tmp_raw_counts = datasets["Maynard_Bivona_2020_NSCLC"].X
datasets["Maynard_Bivona_2020_NSCLC"] = normalize_by_gene_length(
    datasets["Maynard_Bivona_2020_NSCLC"]
)
datasets["Maynard_Bivona_2020_NSCLC"].layers["raw_counts"] = tmp_raw_counts

# %%
datasets["Laughney_Massague_2020_NSCLC"].obs["sex"] = "nan"

# %%
datasets["Maier_Merad_2020_NSCLC"].obs["sex"] = "nan"

# %%
datasets["Madissoon_Meyer_2020_pulmonary-fibrosis"].obs["tissue"] = "lung"
datasets["Madissoon_Meyer_2020_pulmonary-fibrosis"].obs["origin"] = "normal"
datasets["Madissoon_Meyer_2020_pulmonary-fibrosis"].obs["condition"] = "healthy_control"
datasets["Madissoon_Meyer_2020_pulmonary-fibrosis"].obs["sex"] = "nan"
datasets["Madissoon_Meyer_2020_pulmonary-fibrosis"].X.data = np.rint(
    datasets["Madissoon_Meyer_2020_pulmonary-fibrosis"].X.data
)

# %%
datasets["Mayr_Schiller_2020_pulmonary-fibrosis"].obs["sex"] = "nan"
datasets["Mayr_Schiller_2020_pulmonary-fibrosis"] = datasets[
    "Mayr_Schiller_2020_pulmonary-fibrosis"
][
    datasets["Mayr_Schiller_2020_pulmonary-fibrosis"].obs["condition"]
    == "healthy_control",
    :,
]

# %%
datasets["Reyfman_Misharin_2018_pulmonary-fibrosis"].obs["sex"] = "nan"
datasets["Reyfman_Misharin_2018_pulmonary-fibrosis"] = datasets[
    "Reyfman_Misharin_2018_pulmonary-fibrosis"
][
    datasets["Reyfman_Misharin_2018_pulmonary-fibrosis"].obs["condition"]
    == "healthy_control",
    :,
]

# %%
datasets["Travaglini_Krasnow_2020_Lung_SS2"] = datasets[
    "Travaglini_Krasnow_2020_Lung_SS2"
][datasets["Travaglini_Krasnow_2020_Lung_SS2"].obs["tissue"] == "lung", :]
# store true raw counts for later export
tmp_raw_counts = datasets["Travaglini_Krasnow_2020_Lung_SS2"].X
datasets["Travaglini_Krasnow_2020_Lung_SS2"] = normalize_by_gene_length(
    datasets["Travaglini_Krasnow_2020_Lung_SS2"]
)
datasets["Travaglini_Krasnow_2020_Lung_SS2"].layers["raw_counts"] = tmp_raw_counts

# %%
datasets["Zilionis_Klein_2019_NSCLC"] = datasets["Zilionis_Klein_2019_NSCLC"][
    datasets["Zilionis_Klein_2019_NSCLC"].obs["tissue"] == "lung", :
]
datasets["Zilionis_Klein_2019_NSCLC"].obs["sex"] = [
    {"M": "male", "F": "female", "Unknown": "nan"}[s]
    for s in datasets["Zilionis_Klein_2019_NSCLC"].obs["sex"]
]

# %% [markdown]
# ### make patients unique across datasets
#
# Except for the two Travaglini variants - they are the same patients profiled with different platforms

# %%
for dataset_id, adata in datasets.items():
    adata.obs["dataset"] = dataset_id
    adata.obs["patient"] = [
        f'{dataset.replace("_10x", "").replace("_SS2", "")}_{patient}'
        for dataset, patient in zip(adata.obs["dataset"], adata.obs["patient"])
    ]
    datasets[dataset_id] = adata

# %% [markdown]
# ### Validate data

# %%
datasets["Lambrechts_2018_LUAD_6653"].X.data

# %%
for dataset_id, adata in datasets.items():
    print(f"Validating {dataset_id}")
    sanitize_adata(adata)
    validate_adata(adata)

# %% [markdown]
# ## Gene identifier remapping
#
# Use precompiled, static table for gene symbol remapping, since querying MyGene.info requires
# internet connection and is not guaranteed to be reproducible. 

# %%
# datasets_remapped = process_map(remap_gene_symbols, datasets.values(), max_workers=32)
# for dataset_id, dataset in zip(datasets.keys(), datasets_remapped):
#     datasets[dataset_id] = dataset

# %%
# gene_symbol_dict = pd.concat(
#     x.var["original_gene_symbol"]
#     .reset_index()
#     .rename(columns={"index": "gene_symbol", "original_gene_symbol": "alias"})
#     for x in datasets_remapped
# ).drop_duplicates().dropna()
# gene_symbol_dict.to_csv("../../tables/gene_symbol_dict.csv")

# %%
gene_symbol_df = pd.read_csv(
    nxfvars.get("gene_symbol_table", "../../tables/gene_symbol_dict.csv"),
    index_col=False,
)
gene_symbol_dict = {
    alias: symbol
    for alias, symbol in zip(gene_symbol_df["alias"], gene_symbol_df["gene_symbol"])
}

# %%
for dataset_id, tmp_dataset in datasets.items():
    tmp_dataset.var_names = [gene_symbol_dict.get(x, x) for x in tmp_dataset.var_names]

# %% [markdown]
# ### aggregate duplicate gene symbols

# %%
for dataset_id, dataset in datasets.items():
    print(dataset_id)
    datasets[dataset_id] = aggregate_duplicate_gene_symbols(dataset)

# %% [markdown]
# ## Merge data

# %%
obs_all = pd.concat([x.obs for x in datasets.values()], ignore_index=True).reset_index(
    drop=True
)
obs_all = (
    obs_all.loc[
        :,
        MANDATORY_COLS
        + [
            "accession",
            "sampleType",
            "platform",
            "age",
            "tobacco",
            "ethnicity",
            "processing_site",
            "Tissue origins",
            "histology",
            "smoking",
            "pathology",
            "EGFR",
            "tumor_stage",
            "geo_accession",
            "tissue_orig",
            "replicate",
            "race",
            "smoking_status",
            "driver_gene",
            "driver_mutation",
            "secondary_mutation",
            "Notes",
            "stage_at_diagnosis",
            "pathlogy_review",
            "biopsy_date",
            "sort_date",
            "biopsy_type",
            "biopsy_time_status",
            "early_treatment_status",
            "best_response_status",
            "biopsy_timing",
            "analysis",
            "treatment_history",
            "treatment_history_detail",
            "line_of_therapy",
            "treatment_type",
            "treatment",
            "percent_PFS_ref_values",
            "percent.PFS.reference.values",
            "infections",
            "early_bx_day",
            "treatment_start_date",
            "pfs_over_under",
            "pfs_day",
            "pfs_month",
            "date_of_death",
            "stageIII.IV_ca_dx_date",
            "ca_dx_OS",
            "region",
            "location",
            "label",
            "tumor_id",
            "tumor_type",
            "GEO_Sample",
            "biopsy_segment",
            "gsm",
            "characteristics_ch1.7.treatment received prior to surgery (1= treated; 0=untreated)",
        ],
    ]
    .join(dataset_table.set_index("id"), on="dataset")
    .drop_duplicates(ignore_index=False)
    .set_index("sample")
)
# Duplicated doesn't filter out two duplicated rows, don't ask why.
obs_all = obs_all.loc[~obs_all.index.duplicated(), :]

# %%
assert (
    obs_all.index.drop_duplicates().size == obs_all.shape[0]
), "The number of unique samples equals the number of rows"

# %%
merged_all = merge_datasets(datasets.values(), symbol_in_n_datasets=17)

# %%
merged_all.shape

# %% [markdown]
# ### Re-integrate raw counts (not length-scaled) for smartseq2 datasets
#
# This is needed for building a cellxgene-schema-compliant h5ad object 
# in the downstream workflow. 

# %%
# updating values is very slow in sparse matrices.
# We can afford the memory for making this dense temporarily.
merged_all.layers["raw_counts"] = merged_all.X.toarray()

# %%
# the code for "merge all" could have been updated to handle this additional layer.
# however at the point this change was introduced this would have required to patch
# the singularity container again, which would have been more complicated than
# doing this manually here.
for dataset_id in [
    "Guo_Zhang_2018_NSCLC",
    "Maynard_Bivona_2020_NSCLC",
    "Travaglini_Krasnow_2020_Lung_SS2",
]:
    mask = merged_all.obs["dataset"] == dataset_id
    obs_names_merged = merged_all.obs_names[mask]
    obs_names_dataset = obs_names_merged.str.replace("-(\d+)$", "", regex=True)

    common_var = np.intersect1d(merged_all.var_names, datasets[dataset_id].var_names)
    mask_var = merged_all.var_names.isin(common_var)

    # subset dataset in the right order.
    # some cells might not be contained in merged that are in the dataset due to "min batch size"
    tmp_dataset = datasets[dataset_id][obs_names_dataset, common_var]

    assert np.all(
        obs_names_dataset == tmp_dataset.obs_names
    ), "order of indices does not match"
    assert np.all(
        merged_all.var_names[mask_var] == tmp_dataset.var_names
    ), "order of variables does not match"

    merged_all.layers["raw_counts"][np.ix_(mask, mask_var)] = tmp_dataset.layers["raw_counts"].A

# %%
merged_all.layers["raw_counts"] = scipy.sparse.csr_matrix(
    merged_all.layers["raw_counts"]
)

# %% [markdown]
# ## Export all

# %%
merged_all.obs.drop_duplicates().reset_index(drop=True)

# %%
merged_all.write_h5ad(f"{out_dir}/merged_all.h5ad")

# %%
# Some samples drop out due to the min cells threshold. Keep only the remaining samplese in the obs table.
obs_all.loc[merged_all.obs["sample"].unique(), :].to_csv(f"{out_dir}/obs_all.csv")

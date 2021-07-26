# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
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

# %%
out_dir = nxfvars.get("outdir", "/tmp")

# %%
threadpool_limits(int(nxfvars.get("cpus", "8")))

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
dataset_table = pd.read_csv(nxfvars.get("samplesheet", "../../tables/samplesheet_scrnaseq_preprocessing.csv"))

# %%
dataset_table

# %%
dataset_path = nxfvars.get("dataset_path", "../../data/20_qc_norm_scrnaseq/01_qc_and_filtering")
datasets = {
    dataset_id: sc.read_h5ad(
        f"{dataset_id}.qc.h5ad" if dataset_path == "." else f"{dataset_path}/{dataset_id}/{dataset_id}.qc.h5ad"
    )
    for dataset_id in tqdm(dataset_table["id"])
}

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
datasets["Guo_Zhang_2018_NSCLC"] = normalize_by_gene_length(
    datasets["Guo_Zhang_2018_NSCLC"]
)
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
datasets["Maynard_Bivona_2020_NSCLC"] = normalize_by_gene_length(
    datasets["Maynard_Bivona_2020_NSCLC"]
)

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
datasets["Travaglini_Krasnow_2020_Lung_SS2"] = normalize_by_gene_length(
    datasets["Travaglini_Krasnow_2020_Lung_SS2"]
)

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
for dataset_id, adata in datasets.items():
    print(f"Validating {dataset_id}")
    sanitize_adata(adata)
    validate_adata(adata)

# %% [markdown]
# ## Gene identifier remapping

# %%
datasets_remapped = process_map(remap_gene_symbols, datasets.values(), max_workers=32)
for dataset_id, dataset in zip(datasets.keys(), datasets_remapped):
    datasets[dataset_id] = dataset

# %% [markdown]
# ### aggregate duplicate gene symbols
#
# TODO: take care of duplicate genes that were made unique with "var_names_make_unique"

# %%
for dataset_id, dataset in datasets.items():
    print(dataset_id)
    datasets[dataset_id] = aggregate_duplicate_gene_symbols(dataset)

# %% [markdown]
# ## add cell type annotation

# %%
for dataset in datasets.values():
    dataset.obs["cell_type"] = "unknown"

# %% [markdown]
# ## Export all

# %%
obs_all = pd.concat([x.obs for x in datasets.values()]).reset_index(drop=True)

# %%
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
    .drop_duplicates()
    .set_index("sample")
)

# %%
obs_all.to_excel(f"{out_dir}/obs_all.xlsx")

# %%
merged_all = merge_datasets(datasets.values(), symbol_in_n_datasets=17)

# %%
merged_all.shape

# %%
merged_all.obs.drop_duplicates().reset_index(drop=True)

# %%
merged_all.write_h5ad(
    "{out_dir}/merged_all.h5ad"
)

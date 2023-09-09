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
    nxfvars.get("samplesheet", None)
)

# %%
dataset_table

# %%
datasets_annotated = []

# %%
datasets = {
    dataset_id: sc.read_h5ad(
        f"{dataset_id}.qc.h5ad"
    )
    for dataset_id in tqdm(dataset_table["id"])
}

# %%
for dataset_id, adata in tqdm(datasets.items()):
    # Set cell-types of unannotated datasets to "unknown" for scVI
    if "celltype" not in datasets[dataset_id].obs.columns:
        datasets[dataset_id].obs["celltype"] = "unknown"
    
    # Set dataset id
    datasets[dataset_id].obs["dataset"] = dataset_id

    # Validate data
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
    obs_all
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

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
import scanpy as sc
import scanpy_helpers as sh
import numpy as np
import pandas as pd

# %%
adata = sc.read_h5ad("../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad")

# %%
adata_primary = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])
    & (adata.obs["cell_type_major"] != "other"),
    :,
]

# %%
adata_downsampled = sh.deconvolution.balanced_subsample(adata_primary, cell_type_key="cell_type_major", patient_key="patient", n_each=5)

# %%
adata_downsampled

# %%
adata_downsampled.write_h5ad("/home/sturm/Downloads/adata_primary_downsampled.h5ad")

# %%

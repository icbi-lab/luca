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

# %%
import scanpy as sc
from scanpy_helpers.annotation import AnnotationHelper
import warnings
import numpy as np
from nxfvars import nxfvars

warnings.filterwarnings("ignore", category=FutureWarning)

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
input_dir = nxfvars.get(
    "input_dir",
    "../../data/20_integrate_scrnaseq_data/integrate_datasets/27_leiden_umap_nodoublet/",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(f"{input_dir}/adata.doublet_filtered.umap_leiden.h5ad")

# %%
adata.obs["leiden"] = adata.obs["leiden_1.00"]

# %%
sc.pl.umap(adata, color="leiden")

# %%
sc.pl.umap(adata, color="dataset")

# %%
sc.pl.umap(adata, color="dataset", groups=["Maynard_Bivona_2020_NSCLC"], size=8)

# %%
ah.plot_umap(adata, cmap="inferno", size=2)

# %%
sc.pl.umap(adata, color=["GPM6B", "ALB", "HBB"], cmap="inferno", size=2)

# %%
ah.plot_dotplot(adata)

# %%
adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None
)

# %%
sc.pl.umap(
    adata, color="cell_type_predicted", legend_loc="on data", legend_fontoutline=2
)

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "B cell": [8],
    "Epithelial cell": [31, 18, 5, 13, 19, 12, 28, 23, 21, 29],
    "Endothelial cell": [9, 25],
    "Stromal": [14],
    "Granulocytes": [22],
    "Mast cell": [20],
    "Myeloid": [0, 2, 30, 17, 27, 6, 11, 16],
    "Plasma cell": [15],
    "T cell": [7, 1, 4, 3, 10, 24],
    "pDC": [26],
    # they come from brain metastases
    "Neuronal cells": [29],
}

# %%
ah.annotate_cell_types(adata, ct_map)

# %%
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color="dataset", size=1)

# %% [markdown]
# ### Write full adata

# %%
adata.write_h5ad(f"{artifact_dir}/adata_cell_type_coarse.h5ad")

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
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scanpy_helpers.annotation import AnnotationHelper

sc.set_figure_params(figsize=(5, 5))

# %%
adata = sc.read_h5ad(
    "../../data/20_integrate_scrnaseq_data/02_qc_and_filtering/Maynard_Bivona_2020_NSCLC/Maynard_Bivona_2020_NSCLC.qc.h5ad"
)
adata_scvi = sc.read_h5ad(
    "../../data/20_integrate_scrnaseq_data/11_seed_annotations_leiden_umap/integrated_Maynard_Bivona_2020_NSCLC.qc_hvg.umap_leiden.h5ad"
)

assert adata.var_names.is_unique
adata.obsm = adata_scvi.obsm

# %%
adata.obs["leiden"] = adata_scvi.obs["leiden_1.00"]

# %%
adata.X = adata.layers["tpm"]
sc.pp.log1p(adata)

# %%
sc.pl.umap(
    adata, color=["patient", "sample", "CD8A", "CD4", "FOXP3", "leiden"], ncols=3
)

# %% [markdown]
# # Annotation

# %%
ah = AnnotationHelper()

# %%
ah.plot_umap(adata)

# %%
ah.plot_dotplot(adata)

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
adata.obs.columns

# %%
ct_map = {
    "Alevolar cell type 1": [26],
    "Alevolar cell type 2": [15],
    "B cell": [12, 32],
    "Ciliated": [19],
    "Endothelial cell": [8],
    "Epithelial cell": [20, 5, 7, 9, 26, 25, 33, 34],
    "Plasma cell": [6, 39],
    "stromal": [17, 2, 13],
    "Mast cell": [14],
    "Myeloid": [29, 0, 11, 4, 22],
    "NK cell": [18],
    "T cell": [1, 3, 16, 18, 24, 28],
    "pDC": [21],
    "Granulocyte": [10],
    "other C37": [37],
    "other C35": [35],
    "other C38": [38],
    "DC mature": [27],
    "other C36": [36],
    "other C23": [23],
    "other C31": [31],
    "other C30": [30],
}

# %%
ah.annotate_cell_types(adata, ct_map)

# %%
adata.obs["cell_type"].value_counts()

# %% [markdown]
# ## T cell compartment

# %%
adata_t = adata[adata.obs["cell_type"] == "T cell", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_t, leiden_res=0.5)

# %%
ah.plot_umap(adata_t, filter_cell_type=["T cell", "Div"])

# %%
sc.pl.umap(adata_t, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "T cell dividing": [7],
    "T cell CD8": [22, 4, 2, 3, 6],
    "T cell regulatory": [0],
    "T cell CD4": [1, 8, 5],
    "T cell other": [9, 10],
}

# %%
ah.annotate_cell_types(adata_t, cell_type_map=ct_map)

# %%
ah.integrate_back(adata, adata_t)

# %% [markdown]
# ## Myeloid compartment

# %%
adata_m = adata[adata.obs["cell_type"] == "Myeloid", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_m, leiden_res=0.5)

# %%
ah.plot_umap(adata_m, filter_cell_type=["Macro", "Mono", "DC"])

# %%
sc.pl.umap(adata_m, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "cDC": [4, 9],
    "Macrophage": [7, 8, 2, 3, 1, 6],
    "Monocyte": [5, 0],
}

# %%
ah.annotate_cell_types(adata_m, cell_type_map=ct_map)

# %%
ah.integrate_back(adata, adata_m)

# %% [markdown]
# ## Stromal compartment

# %%
adata_s = adata[adata.obs["cell_type"] == "stromal", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_s, leiden_res=0.5)

# %%
ah.plot_umap(adata_s, filter_cell_type=["Fibro", "muscle", "Peri", "Meso"])

# %%
sc.pl.umap(adata_s, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "Smooth muscle cell": [3],
    "Pericyte": [4],
    "Fibroblast": [7, 0,5, 2, 6, 8, 1]
}

# %%
ah.annotate_cell_types(adata_s, ct_map)

# %%
ah.integrate_back(adata, adata_s)

# %% [markdown]
# ## Epithelial compartment

# %%
adata_epi = adata[adata.obs["cell_type"] == "Epithelial cell", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_epi, leiden_res=0.5)

# %%
ah.plot_umap(
    adata_epi,
    filter_cell_type=[
        "Alevolar",
        "Basal",
        "Club",
        "Dividing",
        "Goblet",
        "Ionocyte",
        "Mesothelial",
        "Suprabasal",
    ],
)

# %%
sc.pl.umap(adata_epi, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "Alevolar cell type 1": [8],
    "Alevolar cell type 2": [4],
    "Goblet": [14],
    "Basal": [1, 10, 3, 7, 6, 9],
    "Club": [2,5],
    "Epithelial cell dividing": [11],
    "Epithelial other C12": [12],
    "Mesothelial": [0],
}

# %%
ah.annotate_cell_types(adata_epi, ct_map)

# %%
ah.integrate_back(adata, adata_epi)

# %%
sc.set_figure_params(figsize=(8, 8))
sc.pl.umap(
    adata,
    color="cell_type",
    legend_loc="on data",
    legend_fontoutline=2,
    legend_fontsize=7,
)
sc.set_figure_params(figsize=(5, 5))

# %%
adata.write_h5ad(
    "../../data/30_annotate_scrnaseq_data/maynard_annotated.h5ad", compression="lzf"
)

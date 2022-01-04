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
from nxfvars import nxfvars

sc.set_figure_params(figsize=(5, 5))

# %%
adata = sc.read_h5ad(
    nxfvars.get(
        "adata_qc",
        "../../data/20_integrate_scrnaseq_data/02_qc_and_filtering/Lambrechts_2018_LUAD_6653/Lambrechts_2018_LUAD_6653.qc.h5ad",
    )
)
adata_scvi = sc.read_h5ad(
    nxfvars.get(
        "adata_scvi",
        "../../data/20_integrate_scrnaseq_data/11_seed_annotations_leiden_umap/Lambrechts_2018_LUAD_6653.umap_leiden.h5ad",
    )
)
artifact_dir = nxfvars.get("artifact_dir", "/local/scratch/sturm/")

assert adata.var_names.is_unique
adata.obsm = adata_scvi.obsm

# %%
adata.obs["leiden"] = adata_scvi.obs["leiden_1.00"]

# %%
adata.layers["counts"] = adata.X

# %%
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)

# %%
adata.obs["patient"] = [str(x) for x in adata.obs["patient"]]

# %%
sc.pl.umap(
    adata, color=["patient", "sample", "CD8A", "CD4", "FOXP3", "leiden"], ncols=3
)

# %%
sc.pl.umap(adata, color=["patient", "origin", "condition"])

# %% [markdown]
# ### Annotation

# %%
ah = AnnotationHelper()

# %%
ah.plot_umap(adata)

# %%
ah.plot_dotplot(adata)

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "B cell": [4, 28],
    "Ciliated": [24],
    "Endothelial cell": [13, 29],
    "Epithelial cell": [22, 10, 23, 11, 16, 26],
    "stromal": [21],
    "Myeloid": [0, 30, 25, 6, 3, 8, 14, 12],
    "Mast cell": [17],
    "NK cell": [9],
    "Plasma cell": [15],
    "T cell": [7, 19, 1, 2, 18, 5, 16, 20],
    "pDC": [27],
    "other L31": [31],
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
ah.plot_umap(adata_t, filter_cell_type=["T cell", "Div", "NK"])

# %%
sc.pl.umap(adata_t, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "T cell dividing": [9],
    "NK cell": [8],
    "T cell CD8": [6, 0, 11, 2],
    "T cell regulatory": [3],
    "T cell CD4": [1, 5, 4, 7, 5],
    "other T L10": [10],
}

# %%
ah.annotate_cell_types(adata_t, ct_map)

# %%
ah.integrate_back(adata, adata_t)

# %% [markdown]
# ## myeloid compartment

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
    "DC mature": [10],
    "Monocyte": [6, 7],
    "Macrophage": [11, 0, 2, 8, 4, 1, 5, 9],
    "cDC": [3],
}

# %%
ah.annotate_cell_types(adata_m, ct_map)

# %%
ah.integrate_back(adata, adata_m)

# %% [markdown]
# ## Stromal compartment

# %%
adata_s = adata[adata.obs["cell_type"] == "stromal", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_s)

# %%
ah.plot_umap(adata_s, filter_cell_type=["Fibro", "muscle", "Peri", "Meso"])

# %%
sc.pl.umap(adata_s, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "Mesothelial": [10],
    "Pericyte": [0],
    "Smooth muscle cell": [6],
    "Fibroblast": [1, 3, 8, 4, 7, 9, 5, 2],
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
ah.reprocess_adata_subset_scvi(adata_epi)

# %%
ah.plot_umap(
    adata_epi,
    filter_cell_type=[
        "Alveolar",
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
sc.pl.umap(adata_epi, color=["origin"])

# %%
ct_map = {
    "Alveolar cell type 1": [2],
    "Alveolar cell type 2": [1, 0, 13],
    "Epithelial cell malignant": [3, 10, 5, 8, 9, 14, 6, 7, 12, 4, 11],
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
adata.obs["cell_type"].value_counts().sort_index()

# %%
adata.write_h5ad(f"{artifact_dir}/Lambrechts_2018_LUAD_6653_annotated.h5ad")

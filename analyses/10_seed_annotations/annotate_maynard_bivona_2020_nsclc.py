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
        "../../data/20_integrate_scrnaseq_data/02_qc_and_filtering/Maynard_Bivona_2020_NSCLC/Maynard_Bivona_2020_NSCLC.qc.h5ad",
    )
)
adata_scvi = sc.read_h5ad(
    nxfvars.get(
        "adata_scvi",
        "../../data/20_integrate_scrnaseq_data/11_seed_annotations_leiden_umap/Maynard_Bivona_2020_NSCLC.umap_leiden.h5ad",
    )
)
artifact_dir = nxfvars.get("artifact_dir", "/local/scratch/sturm/")


assert adata.var_names.is_unique
adata.obsm = adata_scvi.obsm

# %%
adata.obs["leiden"] = adata_scvi.obs["leiden_1.00"]

# %%
adata.layers["counts"] = adata.X
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
    "B cell": [12, 32],
    "Ciliated": [19],
    "Endothelial cell": [8],
    "Epithelial cell": [20, 5, 7, 9, 26, 25, 33, 34, 26, 15],
    "Plasma cell": [6, 39],
    "stromal": [17, 2, 13],
    "Mast cell": [14],
    "Myeloid": [29, 0, 11, 4, 22],
    "NK cell": [18],
    "T cell": [1, 3, 16, 18, 24, 28],
    "pDC": [21],
    "Granulocyte": [10],
    "other M37": [37],
    "other M35": [35],
    "other M38": [38],
    "DC mature": [27],
    "other M36": [36],
    "other M23": [23],
    "other M31": [31],
    "other M30": [30],
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
    "T cell dividing": [7],
    "T cell CD8": [22, 4, 2, 3],
    "T cell regulatory": [0],
    "T cell CD4": [1, 8, 5],
    "T cell other": [9, 10],
    "NK cell": [6],
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
    "Fibroblast": [7, 0, 5, 2, 6, 8, 1],
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
sc.pl.umap(adata_epi, color=["origin"])

# %%
sc.pl.umap(adata_epi, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "Alveolar cell type 1": [5, 9],
    "Alveolar cell type 2": [3, 12],
    "Epithelial cell malignant": [14, 13, 0, 2, 11, 8, 4, 6, 1, 7, 10],
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
adata.write_h5ad(f"{artifact_dir}/Maynard_Bivona_2020_NSCLC_annotated.h5ad")

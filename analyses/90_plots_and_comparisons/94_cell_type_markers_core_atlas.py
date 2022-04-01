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
from nxfvars import nxfvars
import altair as alt
from toolz.functoolz import pipe
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy_helpers as sh
import statsmodels.formula.api as smf
from tqdm.contrib.concurrent import process_map
import itertools
import progeny
import dorothea
from natsort import natsorted
from threadpoolctl import threadpool_limits

# %%
alt.data_transformers.disable_max_rows()

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
main_adata = nxfvars.get(
    "main_adata",
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)
epithelial_adata = nxfvars.get(
    "epithelial_adata",
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/epithelial_cells_annotated.h5ad",
)
tumor_adata = nxfvars.get(
    "tumor_adata",
    "../../data/20_build_atlas/annotate_datasets/33_cell_types_epi/artifacts/adata_tumor.h5ad",
)
threadpool_limits(nxfvars.get("cpus", 16))

# %%
adata = sc.read_h5ad(main_adata)

# %%
adata_epi = sc.read_h5ad(epithelial_adata)

# %%
adata_tumor = sc.read_h5ad(tumor_adata)

# %% [markdown]
# # Cell type coarse

# %%
adata.obs["cell_type_coarse"] = pd.Categorical(
    adata.obs["cell_type_coarse"],
    categories=sorted(adata.obs["cell_type_coarse"].unique(), key=lambda x: x.lower()),
)

# %%
sc.pl.dotplot(
    adata,
    groupby="cell_type_coarse",
    var_names={
        "B cell": ["CD79A", "MS4A1"],
        "cDC": ["CLEC9A", "CD1C"],
        "Endothelial cell": ["VWF", "CDH5"],
        "Epithelial cell": ["EPCAM"],
        "Macrophage/Monocyte": ["CD14", "C1QB"],
        "Mast cell": ["TPSB2"],
        "Neutrophils": ["FCGR3B", "CSF3R"],
        "NK cell": ["KLRD1", "GNLY"],
        "pDC": ["CLEC4C"],
        "Plasma cell": ["MZB1"],
        "Stromal": ["COL1A1", "VCAN"],
        "T cell": ["CD3E"],
    },
)

# %% [markdown]
# # cell types fine

# %%
adata.obs["cell_type"] = pd.Categorical(
    adata.obs["cell_type"],
    categories=sorted(adata.obs["cell_type"].unique(), key=lambda x: x.lower()),
)

# %%
marker_dict = {
    "Alveolar cell type 1": ["AGER", "CLDN18"],
    "Alveolar cell type 2": ["SFTPC", "SFTPB", "SFTPA1"],
    "B cell": ["CD19", "CD79A", "MS4A1"],
    "cDC1": ["CLEC9A", "XCR1", "CLNK"],
    "cDC2": ["CD1C", "CLEC10A"],
    "Ciliated": ["PIFO", "FOXJ1", "HYDIN", "CFAP299"],
    "Club": ["SCGB3A1", "SCGB3A2"],
    "DC mature": ["CCR7", "CCL22"],
    "Endothelial cell": ["VWF", "CDH5", "SELE"],
    "Endothelial cell lymphatic": ["CCL21"],
    "Fibroblast": ["PDGFRA", "FAP", "COL1A1"],
    "Fibroblast adventitial": ["MFAP5", "SCARA5"],
    "Fibroblast alveolar": ["ITGA8", "SCN7A"],
    "Macrophage": ["APOE", "CD5L", "MARCO", "C1QB", "TREM2"],
    "Macrophage alveolar": ["FABP4"],
    "Mast cells": ["TPSB2"],
    "Mesothelial": ["MSLN", "CALB2"],
    "Monocyte": ["CD14", "VCAN", "FCN1"],
    "Monocyte conventional": ["S100A12"],
    "Monocyte non-coventional": ["LILRB1", "LILRB2"],
    "Neutrophils": ["FCGR3B", "CSF3R"],
    "NK cell": ["KLRD1", "GNLY"],
    "pDC": ["IL3RA", "CLEC4C"],
    "Pericyte": ["COX4I2", "PDGFRB"],
    "Plasma cell": ["SDC1", "MZB1"],
    "Smooth muscle cell": ["TAGLN", "MYH11"],
    "T cell": ["CD3E"],
    "T cell CD4": ["CD4"],
    "T cell CD8": ["CD8A"],
    "T cell regulatory": ["FOXP3", "IL2RA", "CTLA4"],
    "Dividing": ["MKI67", "CDK1"],
}

# %%
sc.pl.dotplot(adata, groupby="cell_type", var_names=marker_dict)

# %% [markdown]
# # Tumor cells

# %%
# get rid of excluded cells
adata_tumor = adata_tumor[
    adata_tumor.obs_names[adata_tumor.obs_names.isin(adata.obs_names)], :
]
adata_tumor = adata_tumor[
    adata_tumor.obs["cell_type"].str.contains("Tumor cells"), :
].copy()
adata_tumor.obs["cell_type"] = (
    adata_tumor.obs["cell_type"]
    .str.replace("LSCC", "LUSC")
    .str.replace("Tumor cells ", "")
)

# %%
tumor_markers = {
    "Epithelial": ["EPCAM", "B2M"],
    "Alveolar cell type 1": ["AGER", "CLDN18"],
    "Alveolar cell type 2": ["SFTPC", "SFTPB", "SFTPA1"],
    "Ciliated": ["PIFO", "FOXJ1", "HYDIN", "CFAP299"],
    "Club": ["SCGB3A1", "SCGB3A2"],
    "LUAD": ["CD24", "MUC1", "NAPSA", "NKX2-1", "KRT7", "MSLN"],
    "LUSC": ["KRT5", "KRT6A", "TP63", "NTRK2", "SOX2", "KRT17"],
    "NE": ["CHGA", "SYP", "NCAM1", "TUBA1A"],
    "EMT": ["VIM", "SERPINE1", "CDH1"],
    "mitotic": ["MKI67", "TOP2A"],
    "undifferentiated": ["TACSTD2", "AGR2"],
}

# %%
sc.pl.dotplot(adata_tumor, groupby="cell_type", var_names=tumor_markers)

# %%

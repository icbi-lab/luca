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
from nxfvars import nxfvars
import infercnvpy as cnv
from scanpy_helpers.annotation import AnnotationHelper
from scanpy_helpers import de
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from natsort import natsorted
import dorothea
import progeny
# import hierarchical_bootstrapping as hb
from natsort import natsorted
import itertools
from threadpoolctl import threadpool_limits
import numba
from toolz.functoolz import reduce
from operator import or_

# %%
ah = AnnotationHelper()

# %%
sc.set_figure_params(figsize=(4, 4))

# %%
input_adata = nxfvars.get(
    "input_adata",
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/31_cell_types_coarse/by_cell_type/adata_cell_type_coarse_epithelial_cell.umap_leiden.h5ad",
)
threadpool_limits(nxfvars.get("cpus", 1))
numba.set_num_threads(nxfvars.get("cpus", 1))
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
adata = sc.read_h5ad(input_adata)

# %%
adata.obs["leiden"] = adata.obs["leiden_0.50"]

# %%
sc.pl.umap(adata, color="leiden")

# %% [markdown]
# ### Redo UMAP with PAGA

# %%
sc.tl.paga(adata, groups="leiden")

# %%
sc.pl.paga(adata, color="leiden", threshold=0.2)

# %%
sc.tl.umap(adata, init_pos="paga")

# %%
sc.pl.umap(adata, color="dataset")

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    fig = sc.pl.umap(
        adata,
        color=["cell_type", "EPCAM", "dataset", "condition", "origin"],
        cmap="inferno",
        size=2,
        ncols=3,
        return_fig=True,
    )
    fig.savefig(f"{artifact_dir}/overview_epithelial_cluster.pdf", bbox_inches="tight")

# %% [markdown]
# ## Annotation

# %%
sc.pl.umap(adata, color=["GMNC", "FOXJ1", "FOXA3"], cmap="inferno")

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    ah.plot_umap(
        adata,
        filter_cell_type=[
            "Ciliated",
            "Alevolar",
            "Basal",
            "Club",
            "Dividing",
            "Goblet",
            "Ionocyte",
            "Mesothelial",
            "Suprabasal",
        ],
        cmap="inferno",
    )

# %%
ah.plot_dotplot(adata)

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
cell_type_map = {
    "Alevolar cell type 1": [10],
    "Alevolar cell type 2": [1, 14, 8, 12],
    "Ciliated": [6, 15],
    "Club": [3],
    "tumor cells": [2, 9, 13, 17, 11, 4, 0, 16],
}
cell_type_map.update(
    {
        str(k): [k]
        for k in set(map(int, adata.obs["leiden"].values))
        - reduce(or_, map(set, cell_type_map.values()))
    }
)

# %%
ah.annotate_cell_types(adata, cell_type_map)

# %% [markdown]
# ### Subcluster 5

# %%
adata_5 = adata[adata.obs["cell_type"] == "5", :]

# %%
ah.reprocess_adata_subset_scvi(adata_5, leiden_res=0.5)

# %%
sc.pl.umap(adata_5, color=["origin", "dataset"])

# %%
sc.pl.umap(adata_5, color="ROS1")

# %%
ah.plot_umap(adata_5, filter_cell_type=["Alev", "Goblet", "Club"], cmap="inferno")

# %%
ah.plot_dotplot(adata_5)

# %%
sc.pl.umap(adata_5, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%

# %%
ah.annotate_cell_types(
    adata_5,
    {
        "Alevolar cell type 2": [9],
        "Club": [0, 5, 6],
        "ROS1+ healthy epithelial": [2],
        "tumor cells": [8, 4, 3, 1, 7],
    },
)

# %%
ah.integrate_back(adata, adata_5)

# %% [markdown]
# ### Subcluster 7, contains goblet

# %%
adata_7 = adata[adata.obs["cell_type"] == "7", :]

# %%
ah.reprocess_adata_subset_scvi(adata_7, leiden_res=0.5)

# %%
sc.pl.umap(adata_7, color=["origin", "dataset"], wspace=0.5)

# %%
ah.plot_umap(adata_7, filter_cell_type=["Alev", "Goblet", "Club", "Epi"])

# %%
ah.plot_dotplot(adata_7)

# %%
sc.pl.umap(adata_7, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(
    adata_7,
    {"Goblet": [6, 1, 0, 9, 10], "tumor cells": [7, 8, 2, 5, 3, 4]},
)

# %%
ah.integrate_back(adata, adata_7)

# %% [markdown]
# ## Tumor cells

# %%
adata_tumor = adata[adata.obs["cell_type"] == "tumor cells", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_tumor, leiden_res=0.5)

# %%
# if not removed, paga plotting fails
del adata_tumor.uns["leiden_colors"]

# %%
sc.tl.paga(adata_tumor, "leiden")

# %%
sc.pl.paga(adata_tumor, color="leiden", threshold=0.25)

# %%
sc.tl.umap(adata_tumor, init_pos="paga")

# %%
ah.plot_umap(
    adata_tumor, filter_cell_type=["Epi", "Alev", "Goblet", "Club"], cmap="inferno"
)

# %%
print("confounders")
sc.pl.umap(adata_tumor, color=["ALB", "HBB", "HBA1", "GPM6B"])

# %%
print("general")
sc.pl.umap(
    adata_tumor,
    color=["EPCAM", "CDK1", "NEAT1", "MSLN", "origin", "condition", "dataset"],
    cmap="inferno",
    wspace=0.5,
)

print("LUSC")
sc.pl.umap(
    adata_tumor,
    color=["NTRK2", "KRT5", "TP63", "SOX2"],
    cmap="inferno",
)

print("LUAD")
sc.pl.umap(
    adata_tumor,
    color=["MUC1", "NKX2-1", "KRT7", "SFTA2"],
    cmap="inferno",
)

print("EMT")
sc.pl.umap(
    adata_tumor,
    color=["VIM", "NME2", "MIF", "MSLN", "CHGA"],
    cmap="inferno",
)

# %%
sc.pl.umap(adata_tumor, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
adata_tumor_copy = adata_tumor.copy()

# %%
ah.annotate_cell_types(
    adata_tumor,
    cell_type_map={
        "Club": [9],
        "Tumor cells metastasic MSLN+": [14],
        "Tumor cells LSCC mitotic": [3],
        "Tumor cells LSCC": [4, 1],
        "Tumor cells LUAD mitotic": [8],
        "Tumor cells LUAD": [6, 0, 10, 17, 18, 12],
        "Tumor cells EMT": [5],
        "Tumor cells undifferentiated": [7, 11],
        "Tumor cells C2": [2],
        "Tumor cells C13": [13, 15],
        "Tumor cells C16": [16, 19],
    },
)

# %%
ah.integrate_back(adata, adata_tumor)

# %% [markdown]
# ## Write output file

# %%
adata.write_h5ad(f"{artifact_dir}/adata_epithelial.h5ad")
adata_tumor.write_h5ad(f"{artifact_dir}/adata_tumor.h5ad")

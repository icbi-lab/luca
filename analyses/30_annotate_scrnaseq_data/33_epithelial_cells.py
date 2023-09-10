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
from scanpy_helpers.annotation import AnnotationHelper
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from natsort import natsorted

from natsort import natsorted
from threadpoolctl import threadpool_limits
import numba
from toolz.functoolz import reduce
from operator import or_

# %%
ah = AnnotationHelper()

# %%
sc.set_figure_params(figsize=(4, 4))

# %%
# based on Human Lung Cell Atlas
ah2 = AnnotationHelper(
    markers=pd.read_csv(
        nxfvars.get(
            "hlca_markers",
            "../../tables/gene_annotations/hlca_cell_type_signatures.csv",
        )
    )
)

# %%
input_adata = nxfvars.get(
    "input_adata",
    "../../data/20_build_atlas/annotate_datasets/31_cell_types_coarse/by_cell_type/adata_cell_type_coarse_epithelial_cell.umap_leiden.h5ad",
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

# %%
sc.pl.umap(adata, color="dataset")

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    fig = sc.pl.umap(
        adata,
        color=["cell_type", "EPCAM", "dataset", "condition"],
        cmap="inferno",
        size=2,
        ncols=3,
        return_fig=True,
    )
    fig.savefig(f"{artifact_dir}/overview_epithelial_cluster.pdf", bbox_inches="tight")

# %%
ah2.score_cell_types(adata)

# %% [markdown]
# ## Annotation

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    ah.plot_umap(
        adata,
        filter_cell_type=[
            "Ciliated",
            "Alveolar",
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
ah2.plot_umap_scores(
    adata,
    filter_cell_type=[
        "basal",
        "deutero",
        "ciliated",
        "club",
        "goblet",
        "iono",
        "tuft",
        "Neuroend",
        "SMG",
        "AT1",
        "AT2",
    ],
)

# %%
ah.plot_dotplot(adata)

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
cell_type_map = {
    "Alveolar cell type 1": [10],
    "Alveolar cell type 2": [1, 14, 8, 12],
    "Ciliated": [6, 15],
    "transitional club/AT2": [3],
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
sc.pl.umap(adata_5, color=["condition", "dataset"], wspace=0.8)

# %%
sc.pl.umap(adata_5, color=["ROS1"])

# %%
ah.plot_umap(adata_5, filter_cell_type=["Alev", "Goblet", "Club"], cmap="inferno")

# %%
ah.plot_dotplot(adata_5)

# %%
ah2.plot_dotplot_scores(adata_5)

# %%
ah2.plot_umap_scores(adata_5, filter_cell_type=["Club", "Goblet", "AT2", "Cilated"])

# %%
sc.pl.umap(adata_5, color=["n_genes_by_counts", "total_counts"], vmax=[2000, 10000])

# %%
sc.pl.umap(adata_5, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(
    adata_5,
    {
        "Alveolar cell type 2": [9],
        "transitional club/AT2": [0, 5, 6],
        "ROS1+ healthy epithelial": [2],
        "tumor cells": [8, 4, 3, 1, 7],
    },
)

# %%
ah.integrate_back(adata, adata_5)

# %%
sc.tl.rank_genes_groups(adata, groupby="cell_type", groups=["ROS1+ healthy epithelial"], reference="rest")

# %%
sc.pl.rank_genes_groups_dotplot(adata, dendrogram=False, n_genes=50)

# %% [markdown]
# ### Subcluster 7, contains goblet

# %%
adata_7 = adata[adata.obs["cell_type"] == "7", :]

# %%
ah.reprocess_adata_subset_scvi(adata_7, leiden_res=0.5)

# %%
sc.pl.umap(adata_7, color=["dataset"], wspace=0.5)

# %%
ah.plot_umap(adata_7, filter_cell_type=["Alev", "Goblet", "Club", "Epi"])

# %%
ah.plot_dotplot(adata_7)

# %%
ah2.plot_umap_scores(adata_7, filter_cell_type=["Club", "Goblet", "AT2", "Cilated"])

# %%
sc.pl.umap(adata_7, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(
    adata_7,
    {
        "transitional club/AT2": [6,1],
        "Club": [0, 9, 10],
        "tumor cells": [7, 8, 2, 5, 3, 4],
    },
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
sc.pl.umap(adata_tumor, color="leiden_1.00")

# %%
ah.plot_umap(
    adata_tumor,
    filter_cell_type=["Epi", "Alev", "Goblet", "Club", "basal", "tuft", "iono", "Div"],
    cmap="inferno",
)

# %%
print("confounders")
sc.pl.umap(adata_tumor, color=["ALB", "HBB", "HBA1", "GPM6B"])

# %%
markers = {
    "Alveolar type I": ["AGER", "EMP2"],
    "Alveolar type II": ["SFTPA1", "SFTPC"],
    "Club/Goblet": ["SCGB1A1", "SCGB3A1"],
    "Ciliated": ["CAPS", "SNTN"],
    "LUSC": ["KRT5", "KRT6A", "TP63", "NTRK2", "SOX2", "KRT17"],
    "LUAD": ["CD24", "MUC1", "NAPSA", "NKX2-1", "KRT7", "MSLN"],
    "NE": ["CHGA", "SYP", "NCAM1", "TUBA1A"],
    "EMT": ["VIM", "SERPINE1", "CDH1"],
    "Ki67": ["MKI67", "TOP2A"],
    "undifferentiated": ["TACSTD2", "AGR2"],
    "control": ["EPCAM", "B2M"],
}

# %%
sc.pl.umap(
    adata_tumor,
    color=["condition", "dataset"],
    cmap="inferno",
    wspace=0.5,
)

# %%
for name, genes in markers.items():
    print(name)
    sc.pl.umap(
        adata_tumor,
        color=genes,
        cmap="inferno",
        wspace=0.5,
    )

# %%
sc.tl.leiden(adata_tumor, resolution=1, key_added="leiden_1.00")

# %%
sc.pl.dotplot(adata_tumor, groupby="leiden_1.00", var_names=markers)

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.umap(
        adata_tumor,
        color="leiden_1.00",
        legend_loc="on data",
        legend_fontoutline=2,
        size=3,
    )

# %%
adata_tumor_copy = adata_tumor.copy()

# %%
ah.annotate_cell_types(
    adata_tumor,
    cell_type_map={
        "transitional club/AT2": [12],
        "Tumor cells LUAD": [5, 25, 9, 7, 4, 1, 26, 15, 20, 27, 24, 22, 29, 2, 30, 13],
        "Tumor cells LUAD EMT": [6, 18, 23],
        "Tumor cells LUAD MSLN": [21],
        "Tumor cells LUAD NE": [31],
        "Tumor cells NSCLC mixed": [17],
        "Tumor cells LUAD mitotic": [14, 19],
        "Tumor cells LUSC": [0, 8, 28, 16],
        "Tumor cells LUSC mitotic": [3, 11],
        "Hepatocytes": [32],
        "Hemoglobin+": [10],
    },
    column="leiden_1.00",
)

# %%
ah.integrate_back(adata, adata_tumor)

# %% [markdown]
# ## Write output file

# %%
adata.write_h5ad(f"{artifact_dir}/adata_epithelial.h5ad")
adata_tumor.write_h5ad(f"{artifact_dir}/adata_tumor.h5ad")

# %%
# adata_cxg = sc.AnnData(var=adata.raw.var, obs=adata.obs, X=adata.raw.X, obsm=adata.obsm)
# adata_tumor_cxg = sc.AnnData(var=adata_tumor.raw.var, obs=adata_tumor.obs, X=adata_tumor.raw.X, obsm=adata_tumor.obsm)

# %%
# adata_cxg.write_h5ad("/home/sturm/Downloads/2022-03-30_refined_endothelial.h5ad")
# adata_tumor_cxg.write_h5ad("/home/sturm/Downloads/2022-03-30_refined_tumor.h5ad")

# %%

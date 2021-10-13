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
import scanpy as sc
from nxfvars import nxfvars
import infercnvpy as cnv
from scanpy_helpers.annotation import AnnotationHelper

# %%
ah = AnnotationHelper()

# %%
sc.set_figure_params(figsize=(6, 6))

# %%
adata = sc.read_h5ad(
    nxfvars.get(
        "input_adata",
        "../../data/20_integrate_scrnaseq_data/28_annotate_cell_types_coarse_umap/adata_epithelial_cell.umap_leiden.h5ad",
    )
)

# %%
sc.pl.umap(
    adata,
    color=["cell_type", "EPCAM", "dataset", "condition", "origin"],
    cmap="inferno",
    size=2,
    ncols=3,
)

# %%
adata.obs["leiden"] = adata.obs["leiden_0.50"]

# %%
sc.pl.umap(
    adata, color=["leiden_1.00", "leiden_0.50", "leiden_0.75"], legend_loc="on data", palette=""
)

# %%
sc.set_figure_params(figsize=(4,4))

# %%
for c in sorted(adata.obs["leiden"].unique().astype(int)):
    sc.pl.umap(adata, color="leiden", groups=[str(c)], size=1)

# %%
ah.plot_umap(
    adata,
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
    size=1, cmap="inferno"
)

# %%

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
import scvi
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
sc.pl.umap(adata, color="leiden")

# %%
sc.pl.umap(adata, color="dataset")

# %%
sc.pl.umap(adata, color="dataset", groups=["Maynard_Bivona_2020_NSCLC"], size=8)

# %%
ah.plot_umap(adata, cmap="inferno", size=2)

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
    "B cell": [7, 31],
    "Epithelial cell": [18, 21, 11, 13, 3, 12, 30, 25, 27, 29, 32],
    "Endothelial cell": [10, 23],
    "Stromal": [14],
    "Granulocytes": [20],
    "Mast cell": [17],
    "Myeloid": [0, 28, 33, 34, 19, 16, 22, 2, 6, 9, 24],
    "Plasma cell": [15],
    "T cell": [4, 8, 5, 1],
    "pDC": [26],
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

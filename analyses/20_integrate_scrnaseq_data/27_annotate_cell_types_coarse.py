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
    "input_dir", "../../data/20_integrate_scrnaseq_data/26_merge_solo/artifacts"
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(f"{input_dir}/adata.doublet_filtered.h5ad")

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
    "B cell": [7],
    "Epithelial cell": [3, 31, 18, 25, 5, 12, 29, 27, 32, 33, 34, 20],
    "Endothelial cell": [10, 28],
    "Stromal": [16],
    "Granulocytes": [23],
    "Mast cell": [19],
    "Myeloid": [15, 9, 24, 0, 13, 11, 2, 21],
    "Plasma cell": [17],
    "T cell": [6, 26, 22, 1, 8, 4, 14],
    "pDC": [30],
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

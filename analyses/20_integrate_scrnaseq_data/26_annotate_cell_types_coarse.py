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
input_dir = nxfvars.get("input_dir", "/home/sturm/Downloads")
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(f"{input_dir}/adata.doublet_filtered.h5ad")

# %%
sc.pl.umap(adata, color="leiden")

# %%
ah.plot_umap(adata, cmap="inferno", size=2)

# %%
ah.plot_dotplot(adata)

# %%
# adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")
# sc.pp.calculate_qc_metrics(
#     adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None
# )

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "B cell": [7],
    "Ciliated": [19],
    "Endothelial cell": [8],
    "Endothelial cell lymphatic": [30],
    "Epithelial cell": [2, 9, 13, 36, 31, 28, 29, 27, 32, 21, 35, 25, 37],
    "Neutrophil": [24],
    "Mast cell": [18],
    "Myeloid": [1, 16, 33, 6, 5, 12, 10, 17, 38, 26],
    "Stromal": [14],
    "Plasma cell": [15],
    "T cell": [4, 3, 0, 11, 23, 22, 20],
    "pDC": [34],
}

# %%
ah.annotate_cell_types(adata, ct_map)

# %%
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color="dataset", size=1)

# %% [markdown]
# ## Write one adata file per cell-type
# ...  to perform fine-grained annotation in parallel

# %%
for ct in adata.obs["cell_type"].unique():
    ct_filename = ct.replace(" ", "_").lower()
    adata[adata.obs["cell_type"] == ct, :].write_h5ad(f"{artifact_dir}/adata_{ct_filename}.h5ad")

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
input_dir = nxfvars.get("input_dir", "../../data/20_integrate_scrnaseq_data/25_merge_solo/artifacts")
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
# ah.plot_umap(adata, cmap="inferno", size=2)

# %%
ah.plot_dotplot(adata)

# %%
adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None
)

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "B cell": [7,31],
    "Ciliated": [19],
    "Endothelial cell": [10, 33],
    "Endothelial cell lymphatic": [27],
    "Epithelial cell": [8, 29, 26, 28, 3, 20, 30, 12, 35, 24],
    "Neutrophil": [23],
    "Mast cell": [18],
    "Myeloid": [1, 0, 22, 14, 34, 37, 17, 9, 11],
    "Stromal": [15],
    "Plasma cell": [13, 36],
    "T cell": [3, 5, 4, 6, 15, 21, 29, 24],
    "pDC": [32],
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
    tmp_adata = adata[adata.obs["cell_type"] == ct, :].copy()
    tmp_adata.write_h5ad(
        f"{artifact_dir}/adata_{ct_filename}.h5ad"
    )

# %% [markdown]
# ### Write full adata

# %%
adata.write_h5ad(f"{artifact_dir}/adata_cell_type_coarse.h5ad")

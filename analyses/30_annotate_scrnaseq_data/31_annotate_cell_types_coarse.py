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
import warnings
import numpy as np
from nxfvars import nxfvars
import celltypist
from celltypist import models

warnings.filterwarnings("ignore", category=FutureWarning)

# %%
sc.set_figure_params(figsize=(5, 5))


# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(f"adata.doublet_filtered.umap_leiden.h5ad")

# %%
adata_celltypist = adata.copy()  # make a copy of our adata
adata_celltypist.X = adata.layers["raw_counts"]  # set adata.X to raw counts
sc.pp.normalize_per_cell(
    adata_celltypist, counts_per_cell_after=10**4
)  # normalize to 10,000 counts per cell
sc.pp.log1p(adata_celltypist)  # log-transform
# make .X dense instead of sparse, for compatibility with celltypist:
adata_celltypist.X = adata_celltypist.X.toarray()

# %%
models.download_models(model="Cells_Intestinal_Tract.pkl")
model = models.Model.load("Cells_Intestinal_Tract.pkl")

# %% 
predictions = celltypist.annotate(
    adata_celltypist, model=model, majority_voting=True
)
predictions_adata = predictions.to_adata()
predictions_adata

# %%
adata.obs["celltypist_prediction"] = predictions_adata.obs.loc[
    adata.obs.index, "majority_voting"
]

adata.obs["celltypist_conf_score"] = predictions_adata.obs.loc[
    adata.obs.index, "conf_score"
]

# %%
adata.obs["leiden"] = adata.obs["leiden_1.00"]

# %%
sc.pl.umap(adata, color="leiden")

# %%
sc.pl.umap(adata, color="dataset")

# %%
sc.pl.umap(adata, color="celltypist_prediction")

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
adata.obs["cell_type"].value_counts()

# %%
sc.pl.umap(adata, color="dataset", size=1)

# %% [markdown]
# ### Write full adata

# %%
adata.write_h5ad(f"{artifact_dir}/adata_cell_type_coarse.h5ad")

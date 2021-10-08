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
import pandas as pd
from nxfvars import nxfvars
from glob import glob
import scanpy as sc
import numpy as np

sc.settings.set_figure_params(figsize=(10, 10))

# %% [markdown]
# ## Get input data

# %%
input_dir = nxfvars.get("input_dir", "../../data/20_integrate_scrnaseq_data/25_solo/")
adata_path = nxfvars.get(
    "adata_path",
    "../../data/20_integrate_scrnaseq_data/24_scanvi_umap/all.umap_leiden.h5ad",
)
adata_merged_path = nxfvars.get(
    "adata_merged",
    "../../data/20_integrate_scrnaseq_data/21_merge_all/artifacts/merged_all.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
doublet_df = pd.concat(
    [pd.read_csv(x, index_col=0) for x in glob(f"{input_dir}/*.csv")]
)
doublet_df.index = doublet_df.index.str.replace("-0$", "", regex=True)

# %% [markdown]
# Re-integrate all genes (for scVI we were only using HVGs)

# %%
adata_solo = sc.read_h5ad(adata_path)
adata_merged_all = sc.read_h5ad(adata_merged_path)

# %%
assert adata_solo.shape[0] == adata_merged_all.shape[0]

# %%
adata_solo

# %%
adata_merged_all

# %%
adata = adata_merged_all
adata.obs = adata_solo.obs
adata.obsm = adata_solo.obsm
adata.obsp = adata_solo.obsp
adata.uns = adata_solo.uns
adata.var["is_highly_variable"] = ["True" if v in adata_solo.var_names else "False" for v in adata.var_names]

# %% [markdown]
# ## add doublet information to anndata

# %%
adata.obs["doublet_status"] = doublet_df["label"]

# %% [markdown]
# ## doublet statistics
#
# NA = datasets that have been excluded from doublet detection (that is, Smartseq2 datasets)

# %%
adata.obs["doublet_status"].value_counts(dropna=False)

# %% tags=[]
pd.set_option("display.max_rows", 1000)
adata.obs.groupby("sample").apply(
    lambda df: pd.DataFrame().assign(
        doublet_frac=[np.sum(df["doublet_status"] == "doublet") / df.shape[0]],
        n=[df.shape[0]],
    )
).sort_values("doublet_frac")

# %%
sc.pl.umap(adata, color="doublet_status", size=5)

# %% [markdown]
# ## Exclude doublets

# %%
adata_nodoublet = adata[adata.obs["doublet_status"] != "doublet", :]

# %%
sc.pl.umap(adata_nodoublet, color="doublet_status", size=5)

# %% [markdown]
# ## Write output file

# %%
adata.write_h5ad(f"{artifact_dir}/adata.doublet_filtered.h5ad")

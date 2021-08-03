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

# %%
input_dir = nxfvars.get("input_dir", "../../data/20_integrate_scrnaseq_data/23_solo/")
adata_path = nxfvars.get("adata_path", "../../data/20_integrate_scrnaseq_data/22_scvi/integrated_merged_all_all_genes.h5ad")

# %%
doublet_df = pd.concat([pd.read_csv(x, index_col=0) for x in glob(f"{input_dir}/*.csv")])

# %%
doublet_df.index = doublet_df.index.str.replace("-0$", "", regex=True)

# %%
adata = sc.read_h5ad(adata_path)

# %%
adata.obs = adata.obs.join(doublet_df, how="left")

# %%
adata.obs

# %%
np.sum(adata.obs["label"] == "doublet")

# %%
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)

# %%
sc.pl.umap( adata, color="label")

# %%
adata.write_h5ad("/home/sturm/Downloads/adata_doublets.h5ad")

# %%

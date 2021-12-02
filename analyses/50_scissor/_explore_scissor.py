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
import scanpy as sc
import pandas as pd
from nxfvars import nxfvars
sc.settings.set_figure_params(figsize=(5,5))
from pathlib import Path

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
sc.pl.umap(adata, color="cell_type")

# %%
scissor_res_files = Path("/home/sturm/Downloads/scissor/").glob("*_scissor.tsv")

# %%
scissor_ids = [pd.read_csv(x, sep="\t") for x in scissor_res_files]

# %%
scissor_obs = pd.concat(scissor_ids).set_index("cell_id").rename(columns={"Scissor_select": "scissor"})

# %%
adata.obs["scissor"] = scissor_obs["scissor"]

# %%
sc.settings.set_figure_params(figsize=(8, 8))

# %%
sc.pl.umap(adata, color=["scissor", "cell_type"], size=1)

# %%
sc.pl.umap(adata, color=["scissor", "cell_type_coarse"], size=1)

# %%
sc.pl.umap(adata, color=["dataset", "condition", "origin"], size=1)

# %%

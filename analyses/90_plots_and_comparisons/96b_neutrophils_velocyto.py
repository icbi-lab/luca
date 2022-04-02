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
import scvelo as scv
import scanpy as sc
from multiprocessing import Pool
import anndata
import numpy as np
import pandas as pd
from matplotlib import rcParams
from nxfvars import nxfvars
from glob import glob
from tqdm.contrib.concurrent import process_map

rcParams["axes.grid"] = False

# %%
adata_n = sc.read_h5ad(
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/adata_neutrophil_clusters.h5ad"
)
velocyto_dir = nxfvars.get("velocyto_dir", "../../data/11_own_datasets/velocyto/")

# %%
# !ls ../../data/11_own_datasets/velocyto/

# %%
filename_map = {
    "UKIM-V_P1": "Combined_HJNHMDRXX_111723_GCTACGCT_S1_Lanes_final",
    "UKIM-V_P2": "Combined_nsclc242_S1_Lanes_final",
    "UKIM-V_P3": "Combined_124223_final",
    "UKIM-V-2_P4": "Combined_182807_final",
    "UKIM-V-2_P5": "_1_Combined_173026_final",
    "UKIM-V-2_P6": "Combined_182808_final",
    "UKIM-V-2_P7": "Combined_182809_final",
    "UKIM-V-2_P8": "Combined_182810_final",
}


# %%
def _read_scvelo(patient, filename):
    path = list(glob(f"{velocyto_dir}/{filename}_*.loom"))
    assert len(path) == 1
    adata = scv.read_loom(path[0])
    adata.obs["patient"] = patient
    adata.obs["filename"] = filename
    adata.obs_names = [patient + ":" + x.split(":")[1] for x in adata.obs_names]
    return adata


adatas = process_map(_read_scvelo, filename_map.keys(), filename_map.values())

# %%
adata_ukim = adata_n[adata_n.obs["dataset"].str.contains("UKIM"), :].copy()
adata_ukim.obs_names = [f"{patient}:{bc.split('_')[0]}" for patient, bc in zip(adata_ukim.obs["patient"], adata_ukim.obs_names)]

# %%
for ad in adatas:
    ad.var_names_make_unique()

# %%
adata_scvelo = anndata.concat(adatas)

# %%
assert not any(adata_scvelo.obs_names.duplicated())

# %%
# Number of common cells in transcriptomics and scvelo anndata objects.
len(set(adata_scvelo.obs_names)), len(set(adata_ukim.obs_names)), len(
    set(adata_scvelo.obs_names) & set(adata_ukim.obs_names)
)

# %% [markdown]
# ### Preprocess anndata scvelo object

# %%
adata_scvelo = scv.utils.merge(adata_scvelo, adata_ukim)

# %%
adata_ukim.obs["patient"].value_counts()

# %%
adata_scvelo.obs["patient"].value_counts()

# %%
adata_scvelo.shape

# %%
scv.pp.filter_and_normalize(adata_scvelo)

# %%
scv.pp.moments(adata_scvelo, n_pcs=30, n_neighbors=30)

# %%
scv.tl.velocity(adata_scvelo)

# %%
scv.tl.velocity_graph(adata_scvelo)

# %% [markdown]
# ## Plots

# %%
# scv.pl.proportions(adata_scvelo, groupby="patient")

# %%
sc.set_figure_params(figsize=(10, 10))
rcParams["axes.grid"] = False

# %%
ax = sc.pl.embedding(
    adata_ukim,
    basis="umap",
    color="cell_type",
    show=False,
    legend_loc="None",
    size=70,
    alpha=0.4,
)
scv.pl.velocity_embedding_stream(
    adata_scvelo,
    basis="umap",
    color="cell_type",
    #     arrow_color="white",
    legend_loc="right margin",
    ax=ax,
    alpha=0,
    arrow_size=2
)

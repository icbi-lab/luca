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
import pandas as pd
from nxfvars import nxfvars

sc.settings.set_figure_params(figsize=(5, 5))
from pathlib import Path
from scanpy_helpers.annotation import AnnotationHelper
import progeny
import dorothea
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
import altair as alt
import re
import statsmodels.stats
import warnings
import scipy.stats
from tqdm.auto import tqdm
import numpy as np
import statsmodels.formula.api as smf
from tqdm.contrib.concurrent import process_map
import itertools
from operator import or_
from functools import reduce
import scanpy_helpers as sh

# %%
clinical_data_path = nxfvars.get(
    "clincal_data", "../../tables/tcga/clinical_data_for_scissor.tsv"
)

# %%
clinical_data = pd.read_csv(
    clinical_data_path, sep="\t", index_col="TCGA_patient_barcode"
)

# %%
cibersortx_results = pd.read_csv(
    "../../tables/cibersortx/tcga_nsclc/CIBERSORTx_Job2_Results.csv", index_col=0
)
cibersortx_results = cibersortx_results.loc[
    cibersortx_results.index.str.contains("_T$", regex=True)
]

# %%
cibersortx_results.index = cibersortx_results.index.str.replace(
    "(LUAD|LUSC)_", "", regex=True
).str.replace("_(T|N)", "", regex=True)

# %%
ad_immune = sc.AnnData(
    cibersortx_results.drop(
        columns=["P-value", "Correlation", "RMSE", "Absolute score (sig.score)"]
    )
)
ad_immune.obs = clinical_data.loc[ad_immune.obs_names, :]

# %%
ad_immune.obs["patient"] = ad_immune.obs_names

# %% [markdown]
# # Subset anndata objects

# %%
# exclude all non tumor/immune cell-types
# exclude neutrophils, because essentially not present in 10x datasets
EXCLUDE_CELL_TYPES = [
    "Alveolar cell type 1",
    "Alveolar cell type 2",
    "Ciliated",
    "Club",
    "Endothelial cell",
    "Stromal",
    "other",
    "transitional club/AT2",
]

# %%
ad_immune = ad_immune[:, ~ad_immune.var.index.isin(EXCLUDE_CELL_TYPES)]

# %%
ad_immune.var

# %% [markdown]
# ## Infiltration patterns

# %%
# sc.pp.normalize_total(ad_immune, target_sum=1)

# %%
sc.pl.heatmap(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby="type",
    dendrogram=False,
    swap_axes=True,
    cmap="viridis",
    # vmin=-0.25,
    # vmax=0.25
    # # vmin=0,
    # vmax=1,
    standard_scale="var",
)

# %%
sc.pp.neighbors(ad_immune, use_rep="X", n_neighbors=25, metric="correlation")

# %%
sc.pp.scale(ad_immune)

# %%
sc.tl.leiden(ad_immune, resolution=1)

# %%
sc.pl.heatmap(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby="leiden",
    swap_axes=True,
    cmap="bwr",
    vmin=-2.5,
    vmax=2.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %%
ad_immune.obs["immune_type"] = [
    {
        "0": "desert+M",
        "1": "B/CD4/cDC2/alveolar",
        "2": "desert",
        "3": "CD8/cytotoxic",
        "4": "Macrophage",
        "5": "cDC1", 
        "6": "CD8/cytotoxic",
        "7": "plasma",
        "8": "CD8/cytotoxic"
    }[x]
    for x in ad_immune.obs["leiden"]
]

# %%
sc.pl.heatmap(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby="immune_type",
    swap_axes=True,
    cmap="bwr",
    vmin=-2.5,
    vmax=2.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %%
ad_immune.obs.to_csv("/home/sturm/Downloads/stratification_cibersortx.csv")

# %%

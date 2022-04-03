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
from nxfvars import nxfvars
import altair as alt
from toolz.functoolz import pipe
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy_helpers as sh
import statsmodels.formula.api as smf
from tqdm.contrib.concurrent import process_map
import itertools
import progeny
import dorothea
from threadpoolctl import threadpool_limits

# %%
alt.data_transformers.disable_max_rows()

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
main_adata = nxfvars.get(
    "main_adata",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
threadpool_limits(nxfvars.get("cpus", 16))

# %%
adata = sc.read_h5ad(main_adata)

# %% [markdown]
# # Stats and metadata

# %%
adata.X.data.size

# %%
adata.shape[0]

# %%
adata.obs["sample"].nunique()

# %%
adata.obs["platform"].nunique()

# %%
adata.obs["cell_type"].nunique()

# %%
adata.obs["patient"].nunique()

# %%
adata.obs["dataset"].nunique()

# %%
adata.obs["study"].nunique()

# %%
adata.obs.columns

# %%
adata.obs.loc[lambda x: x["origin"].str.contains("tumor")]["patient"].nunique()

# %% [markdown]
# # UMAP plots
#
# UMAP overview plots of the "core" atlas (without projected datasets)

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adata,
        color="cell_type_coarse",
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=1,
        return_fig=True,
        title="",
    )
    # fig.savefig(f"{artifact_dir}/umap_extended_atlas.pdf")

# %%

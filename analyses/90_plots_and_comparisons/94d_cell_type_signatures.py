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

# %% [markdown]
# # Generate signatures
#
# Generate gene signatures for all cell-types that can be used for signature scoring

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
from IPython.display import display_html
import dorothea
from threadpoolctl import threadpool_limits
from matplotlib_venn import venn3
from tqdm.auto import tqdm

# %%
adata_path = nxfvars.get(
    "adata_path",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %% [markdown]
# # Load data

# %%
adata = sc.read_h5ad(adata_path)

# %%
adata.obs["cell_type_major"].unique().tolist()

# %% [markdown]
# # Build signatures

# %%
adata_primary = adata[adata.obs["origin"] == "tumor_primary", :].copy()

# %%
adata_primary_train, adata_primary_test = sh.signatures.train_test_split(
    adata_primary, replicate_col="patient"
)

# %%
pb_primary = sh.pseudobulk.pseudobulk(
    adata_primary,
    groupby=["cell_type_major", "patient", "study"],
)

# %%
sc.pp.normalize_total(pb_primary, target_sum=1e6)
sc.pp.log1p(pb_primary)

# %% [markdown]
# ## Crossvalidation

# %%
# %%time
results = {}
for ct in sorted(adata.obs["cell_type_major"].unique()):
    display_html(f"<h3>working on: {ct}</h3>", raw=True)
    results[ct] = sh.signatures.grid_search_cv(
        adata_primary_train,
        replicate_col="patient",
        label_col="cell_type_major",
        positive_class=ct,
        param_grid={
            "min_fc": list(np.arange(0.5, 3, 0.1)),
            "min_sfc": list(np.arange(0.5, 3, 0.1)),
            "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
        },
    )

# %% [markdown]
# ## Refit and evaluate

# %% [markdown]
# ### Overview

# %%
pd.concat([r.assign(cell_type=ct).iloc[[0]] for ct, r in results.items()]).set_index(
    "cell_type"
).sort_index()

# %%
mcprs = {}
for ct, res in results.items():
    display_html(f"<h3>working on: {ct}</h3>", raw=True)
    mcprs[ct] = sh.signatures.refit_evaluate_plot(
        res, adata_train=adata_primary_train, adata_test=adata_primary_test
    )

# %% [markdown]
# # Evaluate signatures

# %%
sig_table = pd.concat(
    [
        pd.DataFrame().assign(gene_symbol=mcpr.signature_genes).assign(signature=ct)
        for ct, mcpr in mcprs.items()
    ]
).loc[:, ["signature", "gene_symbol"]]

# %%
sig_table.groupby("signature").size().reset_index(name="n")

# %%
assert sig_table["signature"].nunique() == adata.obs["cell_type_major"].nunique()

# %%
for sig, mcpr in tqdm(mcprs.items()):
    genes = mcpr.signature_genes
    sc.tl.score_genes(pb_primary, genes, score_name=sig)
    sc.tl.score_genes(adata, genes, score_name=sig)

# %%
fig = sc.pl.matrixplot(
    pb_primary,
    var_names=sorted(mcprs.keys()),
    cmap="bwr",
    groupby="cell_type_major",
    return_fig=True
)
fig.savefig(f"{artifact_dir}/matrixplot_signature_scores_all_cells.pdf", bbox_inches="tight")

# %%
sc.pl.umap(
    adata,
    color=list(mcprs.keys()),
    cmap="inferno",
    size=1,
    ncols=3,
    frameon=False,
)

# %% [markdown]
# # Export signatures

# %%
sig_table.to_csv(f"{artifact_dir}/cell_type_major_signatures.csv", index=False)

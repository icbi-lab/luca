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
from nxfvars import nxfvars
import altair as alt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import sccoda.util.cell_composition_data as scc_dat
import sccoda.util.comp_ana as scc_ana
import sccoda.util.data_visualization as scc_viz
import scanpy_helpers as sh
import tensorflow as tf

# %%
tf.random.set_seed(0)

# %% [markdown]
# # Input data

# %%
sc.set_figure_params(figsize=(5, 5)) 

# %%
cell_type_column = nxfvars.get("cell_type_column", "cell_type_major")
reference_cell_type = nxfvars.get("reference_cell_type", "Tumor cells")
# Using 500k MCMC iterations. With fewer (tested up to 100k) the results differed
# due to limited random sampling.
mcmc_iterations = nxfvars.get("mcmc_iterations", 10000)
main_adata_file = nxfvars.get(
    "main_adata",
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/full_atlas_neutrophil_clusters.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")


# %%
adata = sc.read_h5ad(main_adata_file)

# %% [markdown]
# # Cell-type fractions

# %%
# only on primary tumor samples;
# exclude datasets with only a single cell-type
frac_by_condition = (
    adata.obs.loc[
        lambda x: (x["origin"] == "tumor_primary")
        & ~x["dataset"].isin(["Guo_Zhang_2018"])
        & x["condition"].isin(["LUAD", "LUSC"])
    ]
    .groupby(["dataset", "condition", "tumor_stage", "patient"])
    .apply(lambda x: x.value_counts(cell_type_column, normalize=False))
    .reset_index(name="n_cells")
    .assign(condition=lambda x: x["condition"].astype(str))
)

# %%
frac_pivot = frac_by_condition.pivot(
    index=["patient", "dataset", "condition", "tumor_stage"],
    columns=cell_type_column,
    values="n_cells",
).reset_index()

# %%
frac_pivot

# %%
data_all = scc_dat.from_pandas(
    frac_pivot, covariate_columns=["patient", "dataset", "condition", "tumor_stage"]
)

# %%
data_all.obs["condition"] = pd.Categorical(data_all.obs["condition"], categories=["LUSC", "LUAD"])

# %%
data_all._sanitize()

# %%
scc_viz.boxplots(data_all, feature_name="condition", figsize=(12, 5))

# %%
scc_viz.stacked_barplot(data_all, feature_name="condition")


# %% [markdown]
# # scCODA

# %%
def run_sccoda(sccoda_data, reference_cell_type, n):
    sccoda_mod = scc_ana.CompositionalAnalysis(
        sccoda_data,
        formula=f"condition + tumor_stage + dataset",
        reference_cell_type=reference_cell_type,
    )
    sccoda_res = sccoda_mod.sample_hmc(num_results=n)
    return sccoda_res


# %%
res_tumor_ref2 = run_sccoda(data_all, reference_cell_type, mcmc_iterations)

# %%
res_tumor_ref2.set_fdr(0.1)

# %%
credible_effects_condition = res_tumor_ref2.credible_effects()[
    "condition[T.LUAD]"
]
credible_effects_stage = res_tumor_ref2.credible_effects()[
    "tumor_stage[T.advanced]"
]

# %%
(
    alt.Chart(
        res_tumor_ref2.effect_df.loc["condition[T.LUAD]"]
        .loc[credible_effects_condition]
        .reset_index(),
        title="condition",
    )
    .mark_bar()
    .encode(
        x=alt.X("Cell Type", sort="y"),
        y="log2-fold change",
        color=alt.Color("Cell Type"),
    )
    | alt.Chart(
        res_tumor_ref2.effect_df.loc["tumor_stage[T.advanced]"]
        .loc[credible_effects_stage]
        .reset_index(),
        title="tumor_stage",
    )
    .mark_bar()
    .encode(
        x=alt.X("Cell Type", sort="y"),
        y="log2-fold change",
        color=alt.Color("Cell Type", scale=sh.colors.altair_scale("cell_type_major")),
    )
).resolve_scale(y="shared", color="shared")

# %%
res_tumor_ref2.effect_df.loc["condition[T.LUAD]"]

# %%
res_tumor_ref2.effect_df.loc["tumor_stage[T.advanced]"]

# %%

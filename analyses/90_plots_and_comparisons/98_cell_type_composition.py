# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:conda-2020-pircher-sccoda]
#     language: python
#     name: conda-env-conda-2020-pircher-sccoda-py
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

# %%
alt.data_transformers.disable_max_rows()

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
main_adata = nxfvars.get(
    "main_adata",
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
)

# %%
adata = sc.read_h5ad(main_adata)

# %%
frac_by_condition = (
    adata.obs.loc[
        lambda x: (x["origin"] == "tumor_primary")
        & ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])
        & x["condition"].isin(["LUAD", "LSCC"])
    ]
    .groupby(["dataset", "condition", "tumor_stage", "patient"])
    .apply(lambda x: x.value_counts("cell_type_major", normalize=False))
    .reset_index(name="n_cells")
    .assign(condition=lambda x: x["condition"].astype(str))
)

# %%
frac_pivot = frac_by_condition.pivot(
    index=["patient", "dataset", "condition", "tumor_stage"],
    columns="cell_type_major",
    values="n_cells",
).reset_index()

# %%
data_all = scc_dat.from_pandas(
    frac_pivot, covariate_columns=["patient", "dataset", "condition", "tumor_stage"]
)

# %%
scc_viz.boxplots(data_all, feature_name="condition", figsize=(12, 5))

# %%
scc_viz.stacked_barplot(data_all, feature_name="condition")


# %%
def run_sccoda(sccoda_data, reference_cell_type, n):
    sccoda_mod = scc_ana.CompositionalAnalysis(
        sccoda_data,
        formula=f"condition + tumor_stage + dataset",
        reference_cell_type=reference_cell_type,
    )
    sccoda_res = sccoda_mod.sample_hmc(num_results=n)
    return sccoda_res


# %% [markdown]
# ## Tumor cells as reference cell-type
#
# Using 500k MCMC iterations. With fewer (tested up to 100k) the results differed 
# due to limited random sampling. 

# %%
res_tumor_ref2 = run_sccoda(data_all, "Tumor cells", 500000)

# %%
credible_effects_condition = res_tumor_ref2.credible_effects(est_fdr=0.1)["condition[T.LUAD]"]
credible_effects_stage = res_tumor_ref2.credible_effects(est_fdr=0.1)["tumor_stage[T.late]"]

# %%
(alt.Chart(
    res_tumor_ref2.effect_df.loc["condition[T.LUAD]"]
    .loc[credible_effects_condition].reset_index(),
    title="condition",
).mark_bar().encode(
    x=alt.X("Cell Type", sort="y"),
    y="log2-fold change",
    color=alt.Color("Cell Type"),
) | alt.Chart(
    res_tumor_ref2.effect_df.loc["tumor_stage[T.late]"]
    .loc[credible_effects_stage]
    .reset_index(),
    title="tumor_stage",
).mark_bar().encode(
    x=alt.X("Cell Type", sort="y"),
    y="log2-fold change",
    color=alt.Color("Cell Type"),
)).resolve_scale(y="shared")

# %%
res_tumor_ref2.effect_df.loc["condition[T.LUAD]"]

# %% [markdown]
# ## Stromal cells as reference
# (for comparison)

# %%
res_tumor_ref2 = run_sccoda(data_all, "Stromal", 500000)

# %%
credible_effects_condition = res_tumor_ref2.credible_effects(est_fdr=0.1)["condition[T.LUAD]"]
credible_effects_stage = res_tumor_ref2.credible_effects(est_fdr=0.1)["tumor_stage[T.late]"]

# %%
(alt.Chart(
    res_tumor_ref2.effect_df.loc["condition[T.LUAD]"]
    .loc[credible_effects_condition].reset_index(),
    title="condition",
).mark_bar().encode(
    x=alt.X("Cell Type", sort="y"),
    y="log2-fold change",
    color=alt.Color("Cell Type"),
) | alt.Chart(
    res_tumor_ref2.effect_df.loc["tumor_stage[T.late]"]
    .loc[credible_effects_stage]
    .reset_index(),
    title="tumor_stage",
).mark_bar().encode(
    x=alt.X("Cell Type", sort="y"),
    y="log2-fold change",
    color=alt.Color("Cell Type"),
)).resolve_scale(y="shared")

# %%
res_tumor_ref2.effect_df.loc["condition[T.LUAD]"]

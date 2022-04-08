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
import scipy.stats
from threadpoolctl import threadpool_limits
import sccoda.util.cell_composition_data as scc_dat
import sccoda.util.comp_ana as scc_ana
import sccoda.util.data_visualization as scc_viz
import tensorflow as tf

# %%
tf.random.set_seed(0)

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
main_adata_file = nxfvars.get(
    "main_adata",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
cell_type_column = nxfvars.get("cell_type_column", "cell_type_major")
reference_cell_type = nxfvars.get("reference_cell_type", "Stromal")
# Using 500k MCMC iterations. With fewer (tested up to 100k) the results differed
# due to limited random sampling.
mcmc_iterations = nxfvars.get("mcmc_iterations", 10000)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
adata = sc.read_h5ad(main_adata_file)

# %% [markdown]
# # Comparison cell-types tumor/normal adjacent

# %%
cell_type_fracs = (
    adata.obs.loc[lambda x: x["origin"].isin(["tumor_primary", "normal_adjacent"])]
    .groupby(["patient", "origin", "study", "dataset", "condition"])
    .apply(
        lambda x: pd.DataFrame()
        .assign(
            n=x[cell_type_column].value_counts(normalize=False),
            frac=x[cell_type_column].value_counts(normalize=True),
        )
        .rename_axis("cell_type")
    )
    .reset_index()
)

# %%
# only keep patient/cell-type combinations that have at least 10 cells in either of the two origins
# only keep patient/cell-type combinations that have both a primary tumor and adjacent normal sample
keep = (
    cell_type_fracs.groupby(["patient", "cell_type"], observed=True)
    .agg(max_cells=("n", np.max), n_origins=("origin", len))
    .loc[lambda x: (x["max_cells"] > 10) & (x["n_origins"] == 2)]
    .reset_index()
    .loc[:, ["patient", "cell_type"]]
)

# %%
plot_cell_type_fracs = (
    cell_type_fracs.merge(keep, how="inner")
    .loc[lambda x: x["cell_type"] != "other"]
    .assign(
        cell_type=lambda x: x["cell_type"]
        .str.replace("Alveolar cell type ", "AT")
        .str.replace("Macrophage", "Macro")
    )
)
for col in ["patient", "origin", "study", "dataset", "cell_type", "condition"]:
    plot_cell_type_fracs[col] = plot_cell_type_fracs[col].astype(str)

# %%
pvalue_df = (
    plot_cell_type_fracs.pivot(
        index=["patient", "cell_type"], columns="origin", values="frac"
    )
    .reset_index()
    .groupby("cell_type")
    .apply(
        lambda x: scipy.stats.ttest_rel(x["normal_adjacent"], x["tumor_primary"]).pvalue
    )
    .reset_index(name="pvalue")
    .pipe(sh.util.fdr_correction)
)

# %%
PROPS = {
    "boxprops": {"facecolor": "none", "edgecolor": "black"},
    "medianprops": {"color": "black"},
    "whiskerprops": {"color": "black"},
    "capprops": {"color": "black"},
}
g = sns.FacetGrid(
    plot_cell_type_fracs,
    col="cell_type",
    aspect=0.6,
    sharey=True,
    legend_out=True,
    col_order=pvalue_df["cell_type"],
    gridspec_kws={"wspace": 0.1, "hspace": 0.2},
)
g.map_dataframe(
    sns.stripplot,
    x="origin",
    y="frac",
    size=5,
    linewidth=1,
    hue="study",
    palette=sh.colors.COLORS.study,
)
g.map_dataframe(sns.boxplot, x="origin", y="frac", color="white", fliersize=0, **PROPS)
g.map_dataframe(
    sns.lineplot,
    x="origin",
    y="frac",
    hue="study",
    ci=None,
    palette=sh.colors.COLORS.study,
)
g.set_xticklabels(rotation=90)
for ax, ct, fdr in zip(g.axes.flatten(), pvalue_df["cell_type"], pvalue_df["fdr"]):
    fdr_str = "FDR<0.01" if fdr < 0.01 else f"FDR={fdr:.2f}"
    ax.set_title(f"{ct}\n{fdr_str}")
g.add_legend()

# %%
PROPS = {
    "boxprops": {"facecolor": "none", "edgecolor": "black"},
    "medianprops": {"color": "black"},
    "whiskerprops": {"color": "black"},
    "capprops": {"color": "black"},
}
g = sns.FacetGrid(
    plot_cell_type_fracs,
    col="cell_type",
    aspect=0.6,
    sharey=True,
    legend_out=True,
    col_order=pvalue_df["cell_type"],
    gridspec_kws={"wspace": 0.1, "hspace": 0.2},
)
g.map_dataframe(
    sns.stripplot,
    x="origin",
    y="frac",
    size=5,
    linewidth=1,
    hue="condition",
    palette=sns.color_palette(),
)
g.map_dataframe(sns.boxplot, x="origin", y="frac", color="white", fliersize=0, **PROPS)
g.map_dataframe(sns.lineplot, x="origin", y="frac", hue="condition", ci=None)
g.set_xticklabels(rotation=90)
for ax, ct, fdr in zip(g.axes.flatten(), pvalue_df["cell_type"], pvalue_df["fdr"]):
    fdr_str = "FDR<0.01" if fdr < 0.01 else f"FDR={fdr:.2f}"
    ax.set_title(f"{ct}\n{fdr_str}")
g.add_legend()

# %%
PROPS = {
    "boxprops": {"facecolor": "none", "edgecolor": "black"},
    "medianprops": {"color": "black"},
    "whiskerprops": {"color": "black"},
    "capprops": {"color": "black"},
}
g = sns.FacetGrid(
    plot_cell_type_fracs.loc[lambda x: x["cell_type"] == "Neutrophils"],
    col="cell_type",
    aspect=0.6,
    sharey=True,
    legend_out=True,
    gridspec_kws={"wspace": 0.1, "hspace": 0.2},
)
g.map_dataframe(
    sns.stripplot,
    x="origin",
    y="frac",
    size=5,
    linewidth=1,
    hue="condition",
    palette=sns.color_palette(),
)
g.map_dataframe(sns.boxplot, x="origin", y="frac", color="white", fliersize=0, **PROPS)
g.map_dataframe(sns.lineplot, x="origin", y="frac", hue="condition", ci=None)
g.set_xticklabels(rotation=90)
g.set_titles("Neutrophils")
g.add_legend()

# %% [markdown]
# # scCODA

# %%
1


# %%
def run_sccoda(sccoda_data, reference_cell_type, n):
    sccoda_mod = scc_ana.CompositionalAnalysis(
        sccoda_data,
        formula=f"origin + patient",
        reference_cell_type=reference_cell_type,
    )
    sccoda_res = sccoda_mod.sample_hmc(num_results=n)
    return sccoda_res


# %%
cell_type_count_pivot = (
    cell_type_fracs.merge(keep, how="inner")
    .pivot(index=["patient", "origin"], columns="cell_type", values="n")
    .fillna(0)
    .reset_index()
    .drop(columns=["other", "Tumor cells"])
)

# %%
cell_type_count_pivot

# %%
sccoda_data = scc_dat.from_pandas(
    cell_type_count_pivot, covariate_columns=["origin", "patient"]
)
sccoda_data._sanitize()

# %%
scc_viz.boxplots(sccoda_data, feature_name="origin", figsize=(12, 5))

# %%
res_tumor_ref2 = run_sccoda(sccoda_data, reference_cell_type, mcmc_iterations)

# %%
res_tumor_ref2.set_fdr(0.1)

# %%
credible_effects_origin = res_tumor_ref2.credible_effects()["origin[T.normal]"]

# %%
alt.Chart(
    res_tumor_ref2.effect_df.loc["origin[T.normal]"]
    .loc[credible_effects_origin]
    .reset_index(),
    title="origin",
).mark_bar().encode(
    x=alt.X("Cell Type", sort="y"),
    y="log2-fold change",
    color=alt.Color("Cell Type"),
)

# %%

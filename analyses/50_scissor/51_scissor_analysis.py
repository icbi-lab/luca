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
import scanpy_helpers as sh

# %%
threadpool_limits(32)

# %%
ah = AnnotationHelper()

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/full_atlas_neutrophil_clusters.h5ad",
)

# %%
scissor_clinical_data = pd.read_csv(
    nxfvars.get("scissor_clinical", "../../tables/tcga/clinical_data_for_scissor.tsv"),
    sep="\t",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
sc.pl.umap(adata, color=["cell_type_coarse", "origin"], wspace=0.8)

# %%
sc.pl.umap(adata, color=["cell_type_neutro", "cell_type_neutro_coarse"], wspace=0.8)

# %% [markdown]
# # Overview clinical data

# %%
scissor_clinical_data.loc[
    :,
    [
        "type",
        "response_to_chemotherapy",
        "tumor_stage",
        "kras_mutation",
        "braf_mutation",
        "egfr_mutation",
        "tp53_mutation",
        "stk11_mutation",
        "stk11_kras_mutation",
    ],
].assign(total=1).groupby("type").agg(sum).astype(int)

# %% [markdown]
# # Load Scissor results

# %%
scissor_res_files = {
    id: Path("../../data/30_downstream_analyses/scissor/scissor_by_sample/").glob(
        f"**/scissor_{id}.tsv"
    )
    for id in [
        "any_tumor_stage",
        "any_response_to_chemotherapy",
        # "any_kras_mutation",
        # "any_braf_mutation",
        # "any_egfr_mutation",
        "any_tumor_type",
        "any_status_time",
        #
        "LUAD_response_to_chemotherapy",
        "LUAD_braf_mutation",
        "LUAD_kras_mutation",
        "LUAD_egfr_mutation",
        "LUAD_tp53_mutation",
        "LUAD_stk11_mutation",
        "LUAD_stk11_kras_mutation",
        "LUAD_tumor_stage",
        "LUAD_random",
        "LUAD_status_time",
        #
        "LUSC_response_to_chemotherapy",
        "LUSC_braf_mutation",
        "LUSC_egfr_mutation",
        "LUSC_tp53_mutation",
        "LUSC_stk11_mutation",
        "LUSC_tumor_stage",
        "LUSC_random",
        "LUSC_status_time",
    ]
}

# %%
scissor_ids = {
    id: [pd.read_csv(x, sep="\t") for x in tmp_files]
    for id, tmp_files in tqdm(scissor_res_files.items())
}

# %%
for id, x in scissor_ids.items():
    assert len(x) > 0, f"{id} has 0 items"

# %%
scissor_obs = {
    f"scissor_{id.lower()}": (
        pd.concat(tmp_ids)
        .set_index("cell_id")
        .rename(columns={"Scissor_select": "scissor"})
    )
    for id, tmp_ids in scissor_ids.items()
}

# %%
for colname, series in scissor_obs.items():
    print(colname)
    adata.obs[colname] = series
    assert adata.obs[colname].value_counts().tolist() == series.value_counts().tolist()

# %%
# pd.set_option("display.max_rows", 100)
# c = "scissor_lusc_status_time"
# adata.obs.loc[lambda x: x["cell_type_coarse"] == "Neutrophils"].groupby(
#     ["patient"]
# ).apply(lambda x: x[c].value_counts(normalize=True, dropna=False)).reset_index().loc[
#     lambda x: ~x["level_1"].isnull()
# ].sort_values(
#     c, ascending=False
# ).head(
#     99
# )

# %%
colname, series = "scissor_lusc_tumor_stage", scissor_obs["scissor_lusc_tumor_stage"]
adata.obs[colname].value_counts().tolist(), series.value_counts().tolist()

# %%
sc.settings.set_figure_params(figsize=(8, 8))

# %%
sc.pl.umap(
    adata,
    color=[
        # "scissor_any_status_time",
        "scissor_luad_status_time",
        "scissor_lusc_status_time",
        "scissor_lusc_tumor_stage",
        "cell_type",
    ],
    size=1,
)

# %% [markdown]
# Scissor+ cells are associated with late stage or with having the corresponding mutation.

# %%
# adata.obs["scissor_tumor_stage"] = [
#     {"scissor+": "late", "scissor-": "early"}.get(x, np.nan)
#     for x in adata.obs["scissor_tumor_stage"]
# ]
# adata.obs["scissor_tumor_type"] = [
#     {"scissor+": "LUSC", "scissor-": "LUAD"}.get(x, np.nan)
#     for x in adata.obs["scissor_tumor_type"]
# ]
# adata.obs["scissor_kras_mutation"] = [
#     {"scissor+": "KRAS mut", "scissor-": "no KRAS mut"}.get(x, None)
#     for x in adata.obs["scissor_kras_mutation"]
# ]
# adata.obs["scissor_braf_mutation"] = [
#     {"scissor+": "BRAF mut", "scissor-": "no BRAF mut"}.get(x, None)
#     for x in adata.obs["scissor_braf_mutation"]
# ]
# adata.obs["scissor_egfr_mutation"] = [
#     {"scissor+": "EGFR mut", "scissor-": "no EGFR mut"}.get(x, None)
#     for x in adata.obs["scissor_egfr_mutation"]
# ]

# %%
adata_primary = adata[
    (adata.obs["origin"] == "tumor_primary"),
    # & ~adata.obs["dataset"].str.startswith("UKIM"),
    :,
].copy()

# %%
sc.pl.umap(adata_primary, color=["dataset", "condition", "origin"], size=1)

# %%
scissor_cols = adata_primary.obs.columns[
    adata_primary.obs.columns.str.startswith("scissor_")
]

# %%
for var in scissor_cols:
    with plt.rc_context({"figure.figsize": (6, 6), "figure.dpi": 300}):
        sc.pl.umap(
            adata_primary,
            color=var,
            size=2,
            palette=["#ca0020", "#0571b0"][::-1],
            frameon=False,
        )


# %%
def _scissor_test(df):
    c1, c2 = [x for x in df.columns if x != "nan"]
    _, p = scipy.stats.wilcoxon(df[c1].values, df[c2].values, zero_method="zsplit")
    return pd.Series(
        {
            c1: df[c1].mean(),
            c2: df[c2].mean(),
            "pvalue": p,
            "log2_ratio": np.log2(df[c1].mean()) - np.log2(df[c2].mean()),
        }
    )


def scissor_by_group(
    adata,
    *,
    groupby=["cell_type_major", "patient"],
    scissor_col,
    cell_cutoff=10,
):
    """Aggregate scissor scores first by patient, then by a grouping variable"""
    obs = adata.obs.copy()
    # convert to str that nans are counted
    obs[scissor_col] = obs[scissor_col].astype(str)
    df_grouped = (
        obs.groupby(groupby, observed=True)[scissor_col]
        .value_counts(normalize=False)
        .reset_index(name="n")
    ).pivot_table(
        values="n",
        columns=scissor_col,
        index=groupby,
        fill_value=0,
        aggfunc=np.sum,
    )
    # only consider samples with at least 10 cells of each type
    df_grouped = df_grouped.loc[
        lambda x: x.apply(lambda row: np.sum(row), axis=1) > cell_cutoff
    ]
    df_grouped += 1  # add pseudocount
    df_grouped = df_grouped.apply(lambda row: row / np.sum(row), axis=1)  # normalize

    df_grouped = (
        df_grouped.groupby(groupby[0]).apply(_scissor_test).pipe(sh.util.fdr_correction)
    )

    return df_grouped


# %%
def plot_scissor_df(df, *, title="scissor"):
    """Plot the result of scissor_by_group as a bar chart"""
    up, down = [x for x in df.columns[1:] if x != "nan"]
    order = df.sort_values(up)["cell_type_major"].values.tolist()
    return alt.vconcat(
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("cell_type_major", sort=order, axis=None),
            y=alt.Y(up, scale=alt.Scale(domain=[0, 1])),
            # color=alt.Color(
            #     "gini_better", scale=alt.Scale(scheme="magma", reverse=False)
            # ),
        )
        .properties(height=100, title=title),
        alt.Chart(df.assign(**{down: lambda x: -x[down]}))
        .mark_bar()
        .encode(
            x=alt.X("cell_type_major", sort=order),
            y=alt.Y(down, scale=alt.Scale(domain=[-1, 0]))
            # color=alt.Color(
            #     "gini_worse", scale=alt.Scale(scheme="magma", reverse=False)
            # ),
        )
        .properties(height=100),
        spacing=0,
    )


# %%
def plot_scissor_df_ratio(
    df, *, title="scissor", fdr_cutoff=0.01, groupby="cell_type_major"
):
    """Plot the result of scissor_by_group as a bar chart"""
    df = df.loc[lambda x: x["fdr"] < fdr_cutoff].copy().reset_index(drop=False)
    # print(df)
    up, down = df.columns[:2]
    # df["log2_ratio"] = np.log2(df[up]) - np.log2(df[down])
    order = df.sort_values("log2_ratio")[groupby].values.tolist()
    max_ = np.max(np.abs(df["log2_ratio"]))
    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X(groupby, sort=order),
            y=alt.Y("log2_ratio", scale=alt.Scale(domain=[-5, 5])),
            color=alt.Color(
                "log2_ratio", scale=alt.Scale(scheme="redblue", domain=[-max_, max_])
            )
            # color=alt.Color(
            #     "gini_better", scale=alt.Scale(scheme="magma", reverse=False)
            # ),
        )
        .properties(height=100, title=title)
    )


# %%
scissor_dfs = {
    k: scissor_by_group(adata_primary, scissor_col=k, cell_cutoff=1)
    for k in scissor_cols
}

# %%
for col, df in scissor_dfs.items():
    plot_scissor_df_ratio(df, title=col).display()

# %%
scissor_dfs = {
    k: scissor_by_group(
        adata_primary,
        scissor_col=k,
        groupby=["cell_type_neutro_coarse", "patient"],
        cell_cutoff=1,
    )
    for k in scissor_cols
}

# %%
for col, df in scissor_dfs.items():
    plot_scissor_df_ratio(df, title=col, groupby="cell_type_neutro_coarse").display()

# %%
scissor_dfs = {
    k: scissor_by_group(
        adata_primary,
        scissor_col=k,
        groupby=["cell_type_neutro", "patient"],
        cell_cutoff=1,
    )
    for k in scissor_cols
}

# %%
for col, df in scissor_dfs.items():
    plot_scissor_df_ratio(df, title=col, groupby="cell_type_neutro").display()

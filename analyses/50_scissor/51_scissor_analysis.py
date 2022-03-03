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
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
sc.pl.umap(adata, color=["cell_type_coarse", "origin"])

# %%
scissor_res_files = {
    id: Path("../../data/30_downstream_analyses/scissor/scissor_by_sample/").glob(
        f"**/scissor_{id}.tsv"
    )
    for id in [
        # "any_status_time",
        "any_tumor_type",
        # "any_tumor_stage",
        # "any_kras_mutation",
        # "any_braf_mutation",
        # "any_egfr_mutation",
        "LUAD_kras_mutation",
        "LUAD_braf_mutation",
        "LUAD_egfr_mutation",
        "LUAD_tp53_mutation",
        "LUAD_tumor_stage",
        "LUAD_random",
        "LUAD_status_time",
        "LUSC_kras_mutation",
        "LUSC_braf_mutation",
        "LUSC_egfr_mutation",
        "LUSC_tumor_stage",
        "LUSC_status_time",
        "LUSC_tp53_mutation",
        "LUSC_random",
    ]
}

# %%
scissor_ids = {
    id: [pd.read_csv(x, sep="\t") for x in tmp_files]
    for id, tmp_files in scissor_res_files.items()
}

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
#     {"scissor+": "LSCC", "scissor-": "LUAD"}.get(x, np.nan)
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
    (adata.obs["origin"] == "tumor_primary")
    & ~adata.obs["dataset"].str.startswith("UKIM"),
    :,
].copy()

# %%
sc.pl.umap(adata_primary, color=["dataset", "condition", "origin"], size=1)

# %%
for var in [
    "scissor_any_tumor_type",
    "scissor_luad_kras_mutation",
    "scissor_luad_braf_mutation",
    "scissor_luad_egfr_mutation",
    "scissor_luad_tp53_mutation",
    "scissor_luad_tumor_stage",
    "scissor_luad_random",
    "scissor_luad_status_time",
    "scissor_lusc_kras_mutation",
    "scissor_lusc_braf_mutation",
    "scissor_lusc_egfr_mutation",
    "scissor_lusc_tumor_stage",
    "scissor_lusc_status_time",
    "scissor_lusc_tp53_mutation",
    "scissor_lusc_random",
]:
    with plt.rc_context({"figure.figsize": (6, 6), "figure.dpi": 120}):
        sc.pl.umap(
            adata_primary,
            color=var,
            size=2,
            palette=["#ca0020", "#0571b0"][::-1],
            frameon=False,
        )

# %%
# def conf_int(data, confidence=0.95):
#     a = 1.0 * np.array(data)
#     n = len(a)
#     m, se = np.mean(a), scipy.stats.sem(a)
#     h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
#     return h


def _scissor_test(df):
    c1, c2 = [x for x in df.columns if x != "nan"]
    _, p = scipy.stats.ttest_rel(df[c1].values, df[c2].values)
    return pd.Series(
        {
            c1: df[c1].mean(),
            c2: df[c2].mean(),
            "pvalue": p,
            "log2_ratio": np.log2(df[c1].mean()) - np.log2(df[c2].mean()),
        }
    )


def scissor_by_group(
    adata, *, groupby=["cell_type_major", "patient"], scissor_col, adatas_for_gini=None
):
    """Aggregate scissor scores first by patient, then by a grouping variable

    Parameters
    ----------
    adatas_for_gini
        Set this to a dict {groupby: adata} with an AnnData object
        for each group with a `leiden` clustering. This will
        be used to obtain gini index.
    """
    obs = adata_primary.obs.copy()
    # convert to str that nans are counted
    obs[scissor_col] = obs[scissor_col].astype(str)
    df_grouped = (
        (
            obs.groupby(groupby, observed=True)[scissor_col]
            .value_counts(normalize=True)
            .reset_index(name="frac")
        )
        .pivot_table(
            values="frac",
            columns=scissor_col,
            index=["cell_type_major", "patient"],
            fill_value=0,
            aggfunc=np.mean,
        )
        .groupby("cell_type_major")
        .apply(_scissor_test)
    ).pipe(sh.util.fdr_correction)

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
def plot_scissor_df_ratio(df, *, title="scissor", fdr_cutoff=0.01):
    """Plot the result of scissor_by_group as a bar chart"""
    df = df.loc[lambda x: x["fdr"] < fdr_cutoff].copy().reset_index(drop=False)
    # print(df)
    up, down = df.columns[:2]
    # df["log2_ratio"] = np.log2(df[up]) - np.log2(df[down])
    order = df.sort_values("log2_ratio")["cell_type_major"].values.tolist()
    max_ = np.max(np.abs(df["log2_ratio"]))
    return (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("cell_type_major", sort=order),
            y=alt.Y("log2_ratio", scale=alt.Scale(domain=[-13, 13])),
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
    k: scissor_by_group(adata, scissor_col=k)
    for k in [
        "scissor_any_tumor_type",
        "scissor_luad_kras_mutation",
        "scissor_luad_braf_mutation",
        "scissor_luad_egfr_mutation",
        "scissor_luad_tp53_mutation",
        "scissor_luad_tumor_stage",
        "scissor_luad_random",
        "scissor_luad_status_time",
        "scissor_lusc_kras_mutation",
        "scissor_lusc_braf_mutation",
        "scissor_lusc_egfr_mutation",
        "scissor_lusc_tumor_stage",
        "scissor_lusc_status_time",
        "scissor_lusc_tp53_mutation",
        "scissor_lusc_random",
    ]
}

# %%
for col, df in scissor_dfs.items():
    plot_scissor_df_ratio(df, title=col).display()


# %% [markdown]
# ---
# ### gini per group

# %%
def subcluster_adata(tmp_adata):
    sc.pp.neighbors(tmp_adata, use_rep="X_scANVI")
    sc.tl.leiden(tmp_adata, resolution=0.5)
    return adata


# %%
adatas = []
cell_types = adata_primary.obs["cell_type"].unique()
for cell_type in cell_types:
    adatas.append(adata_primary[adata_primary.obs["cell_type"] == cell_type, :].copy())

# %%
adatas = process_map(subcluster_adata, adatas)


# %%
def scissor_gini(tmp_adata):
    # tmp_adata = adata[adata.obs["cell_type"] == cell_type, :].copy()
    sc.pp.neighbors(tmp_adata, use_rep="X_scANVI")
    sc.tl.leiden(tmp_adata, resolution=0.5)
    fractions = (
        tmp_adata.obs.groupby("leiden")["scissor"]
        .value_counts(normalize=True)
        .reset_index(name="frac")
        .pivot_table(values="frac", columns="scissor", index="leiden", fill_value=0)
        .reset_index()
    )
    try:
        gini_better = gini_index(fractions["better survival"].values)
    except KeyError:
        gini_better = 0
    try:
        gini_worse = gini_index(fractions["worse survival"].values)
    except KeyError:
        gini_worse = 0
    return gini_better, gini_worse


# %%

# %%
gini_better, gini_worse = zip(*gini_res)

# %%
gini_better = pd.Series(gini_better, index=cell_types)
gini_worse = pd.Series(gini_worse, index=cell_types)

# %%
scissor_per_cell_type.set_index("cell_type", inplace=True)

# %%
scissor_per_cell_type["gini_better"] = gini_better
scissor_per_cell_type["gini_worse"] = gini_worse
scissor_per_cell_type = scissor_per_cell_type.reset_index()

# %% [markdown]
# ---

# %% [markdown]
# # Variability within cell-types

# %%
# %%time
hb.tl.pseudobulk(progeny_epi, groupby=["leiden", "patient"], aggr_fun=np.mean)

# %%
progeny_epi_pb_leiden = hb.tl.pseudobulk(
    progeny_epi, groupby=["patient", "leiden"], aggr_fun=np.mean
)

# %%
progeny_epi_pb = hb.tl.pseudobulk(progeny_epi, groupby=["patient"], aggr_fun=np.mean)

# %%
progeny_epi_pb.X -= np.min(progeny_epi_pb.X)
progeny_epi_pb_leiden.X -= np.min(progeny_epi_pb_leiden.X)

# %%
res = []
patients = progeny_epi_pb_leiden.obs["patient"].unique()
for patient in tqdm(patients):
    tmp_x = progeny_epi_pb_leiden.X[progeny_epi_pb_leiden.obs["patient"] == patient, :]
    res.append(np.apply_along_axis(gini_index, 0, tmp_x))

# %%
gini_within = np.mean(np.vstack(res), axis=0)

# %%
gini_between = np.apply_along_axis(gini_index, 0, progeny_epi_pb.X)

# %%
df_to_plot = (
    pd.DataFrame(gini_within.T, index=progeny_epi.var_names, columns=["gini_within"])
    .join(
        pd.DataFrame(
            gini_between, index=progeny_epi.var_names, columns=["gini_between"]
        )
    )
    .reset_index()
)

# %%
df_to_plot

# %%
import altair as alt

# %%
alt.Chart(df_to_plot).mark_point().encode(
    y="gini_within", x="gini_between", color=df_to_plot.columns[0], tooltip="index"
)

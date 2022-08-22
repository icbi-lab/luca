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
# %config InlineBackend.figure_formats = ['png']
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
import threadpoolctl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st
import matplotlib
import scanpy_helpers as sh
import altair as alt
import numpy as np
from nxfvars import nxfvars

# %%
threadpoolctl.threadpool_limits(20)

# %%
matplotlib.rcParams.update({"font.size": 16})
matplotlib.rcParams["figure.dpi"] = 300

# %%
adata = sc.read_h5ad(
    nxfvars.get(
        "main_adata",
        "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
    )
)

# %% [markdown]
# ## Total cell-type fractions

# %%
df = (
    adata.obs.loc[lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])][
        "cell_type_coarse"
    ]
    .value_counts(normalize=True)
    .reset_index(name="fraction")
    .rename(columns={"index": "cell_type"})
)
alt.Chart(df).mark_bar().encode(
    color=alt.Color(
        "cell_type", scale=sh.colors.altair_scale("cell_type_coarse"), legend=None
    ),
    x="fraction",
    y=alt.Y("cell_type", sort="-x"),
)

# %%
df

# %% [markdown]
# ## Neutrophil fractions per platform (of all cells)

# %%
df = (
    adata.obs.loc[lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])]
    .groupby(
        [
            "platform",
            "dataset",
        ]
    )["cell_type_coarse"]
    .value_counts(normalize=True)
    .reset_index(name="fraction")
    .loc[lambda x: x["level_2"] == "Neutrophils"]
)

# %%
order = (
    df.groupby("platform")
    .agg(np.mean)
    .sort_values("fraction", ascending=False)
    .index.tolist()
)
alt.Chart(df).mark_bar().encode(
    color=alt.Color("platform", scale=sh.colors.altair_scale("platform"), legend=None),
    x=alt.X("mean(fraction)", title="mean neutrophil fraction across all datasets"),
    y=alt.Y("platform", sort=order),
)
# + alt.Chart(df).mark_errorbar(extent="ci").encode(
#     x=alt.X("fraction"),
#     y=alt.Y("platform", sort=order),
# )

# %% [markdown]
# ## Neutrophil fractions per platform (of total Leukocytes)

# %%
df = (
    adata.obs.loc[lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])]
    .loc[
        lambda x: ~x["cell_type_coarse"].isin(
            ["Epithelial cell", "Endothelial cell", "Stromal"]
        )
    ]
    .groupby(
        [
            "platform",
            "dataset",
        ]
    )["cell_type_coarse"]
    .value_counts(normalize=True)
    .reset_index(name="fraction")
    .loc[lambda x: x["level_2"] == "Neutrophils"]
)

# %%
order = (
    df.groupby("platform")
    .agg(np.mean)
    .sort_values("fraction", ascending=False)
    .index.tolist()
)
alt.Chart(df).mark_bar().encode(
    color=alt.Color("platform", scale=sh.colors.altair_scale("platform"), legend=None),
    x=alt.X("mean(fraction)", title="mean neutrophil fraction across all datasets"),
    y=alt.Y("platform", sort=order),
)
# + alt.Chart(df).mark_errorbar(extent="ci").encode(
#     x=alt.X("fraction"),
#     y=alt.Y("platform", sort=order),
# )

# %% [markdown]
# ## Platform fractions per cell-type

# %%
df = (
    (
        adata.obs.loc[
            lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])
        ]["platform"].value_counts(normalize=True)
    )
    .reset_index(name="fraction")
    .rename(columns={"index": "platform"})
    .assign(y="all cells")
)

c_vline = (
    alt.Chart(df.loc[lambda x: x["platform"] == "10x"])
    .mark_rule(color="black", strokeWidth=1, strokeDash=[3, 3])
    .encode(x="fraction")
)


c1 = (
    alt.Chart(df)
    .mark_bar()
    .encode(
        x=alt.X(
            "fraction",
            axis=None,
        ),
        order=alt.Order("fraction", sort="descending"),
        color=alt.Color(
            "platform",
            scale=sh.colors.altair_scale("platform"),
            sort=alt.EncodingSortField("fraction"),
        ),
        y=alt.Y("y", title=None),
    )
) + c_vline

df = (
    adata.obs.loc[lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])]
    .groupby("cell_type_coarse")["platform"]
    .value_counts(normalize=True)
    .reset_index(name="fraction")
    .rename(columns={"level_1": "platform"})
    .set_index("platform")
    .assign(order=df.set_index("platform")["fraction"])
    .reset_index(drop=False)
)
c2 = (
    alt.Chart(df)
    .mark_bar()
    .encode(
        y=alt.Y(
            "cell_type_coarse",
            title=None,
            sort=df.loc[lambda x: x["platform"] == "10x"]
            .sort_values("fraction", ascending=False)["cell_type_coarse"]
            .tolist(),
        ),
        x=alt.X(
            "fraction",
            sort=df.sort_values("fraction")["platform"].tolist(),
            scale=alt.Scale(domain=[0, 1]),
        ),
        color=alt.Color(
            "platform",
            scale=sh.colors.altair_scale("platform"),
            sort=alt.EncodingSortField("order"),
        ),
        order=alt.Order("order", sort="descending"),
    )
) + c_vline

c1 & c2

# %%
df.loc[lambda x: x["cell_type_coarse"] == "Neutrophils"]

# %%
ukimv_fracs = adata.obs.loc[lambda x: x["study"] == "UKIM-V", "cell_type_coarse"].value_counts(normalize=True)

# %%
ukimv_fracs

# %%
ukimv_fracs["Neutrophils"]

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true tags=[]
# # mRNA content
#
# Need to compute ratios, because the baseline difference between datasets and platforms is very high.

# %%
rel_counts = (
    adata.obs.groupby(["dataset", "cell_type_coarse"])
    .agg(total_counts=("total_counts", "median"))
    .reset_index()
    .groupby("dataset")
    .apply(
        lambda x: x.assign(
            rel_counts=np.log2(x["total_counts"])
            - np.log2(
                x.loc[x["cell_type_coarse"] == "Epithelial cell", "total_counts"].values
            )
        )
    )
)

# %%
rel_counts

# %%
order = (
    rel_counts.groupby("cell_type_coarse")
    .mean()
    .sort_values("rel_counts")
    .index.tolist()
)
(
    alt.Chart(
        rel_counts,
        title="Mean detected counts per cell-type, relative to Epithelial cells",
    )
    .mark_bar()
    .encode()
    .encode(
        x=alt.X("cell_type_coarse", sort=order),
        y=alt.Y("mean(rel_counts):Q", title="log2(ratio)"),
        color=alt.Color(
            "cell_type_coarse",
            scale=sh.colors.altair_scale("cell_type_coarse"),
            legend=None,
        ),
    )
    + alt.Chart(rel_counts)
    .mark_errorbar(extent="ci")
    .encode(
        x=alt.X("cell_type_coarse", sort=order),
        y=alt.Y("rel_counts", title="log2(ratio)"),
    )
)

# %%
adata.obs.columns

# %%
for title, tmp_adata in {
    "counts per platform": adata,
    "counts per platform (Epithelial cells)": adata[
        adata.obs["cell_type_coarse"] == "Epithelial cell", :
    ],
}.items():
    counts_per_platform = (
        tmp_adata.obs.groupby(["sample", "platform"], observed=True)
        .agg(total_counts=("total_counts", "median"))
        .reset_index()
    )
    order = (
        counts_per_platform.groupby("platform")
        .median()
        .sort_values("total_counts")
        .index.tolist()
    )
    alt.Chart(counts_per_platform, title=title).mark_boxplot().encode(
        x=alt.X("platform", sort=order[::-1]),
        y=alt.Y("total_counts", scale=alt.Scale(type="log")),
        color=alt.Color(
            "platform",
            scale=sh.colors.altair_scale("platform"),
            legend=None,
        ),
    ).display()

# %%
adata.obs.loc[:, ["platform_fine", "patient"]].drop_duplicates().groupby(
    "platform_fine"
).size()

# %%
for title, tmp_adata in {
    "counts per platform": adata,
    "counts per platform (Epithelial cells)": adata[
        adata.obs["cell_type_coarse"] == "Epithelial cell", :
    ],
}.items():
    counts_per_platform = (
        tmp_adata.obs.groupby(["sample", "platform_fine"], observed=True)
        .agg(total_counts=("total_counts", "median"))
        .reset_index()
    )
    order = (
        counts_per_platform.groupby("platform_fine")
        .median()
        .sort_values("total_counts")
        .index.tolist()
    )
    alt.Chart(counts_per_platform, title=title).mark_boxplot().encode(
        x=alt.X("platform_fine", sort=order[::-1]),
        y=alt.Y("total_counts", scale=alt.Scale(type="log")),
        color=alt.Color(
            "platform_fine",
            # scale=sh.colors.altair_scale("platform"),
            legend=None,
        ),
    ).display()

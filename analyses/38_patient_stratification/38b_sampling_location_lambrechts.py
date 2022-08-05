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
ah = AnnotationHelper()

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(path_adata)

# %%
sc.pl.umap(adata, color=["cell_type_major", "origin"], ncols=1)

# %% [markdown]
# # Subset anndata objects

# %%
adata.obs["cell_type_major"].value_counts().sort_index()

# %%
# exclude all non tumor/immune cell-types
# exclude neutrophils, because essentially not present in 10x datasets
EXCLUDE_CELL_TYPES = [
    "Alveolar cell type 1",
    "Alveolar cell type 2",
    "Neutrophils",
    "Ciliated",
    "Club",
    "Endothelial cell",
    "Stromal",
    "other",
    "transitional club/AT2",
]

# %%
# For the patient stratification, treat the two batches of the UKIM-V dataset as one
adata.obs["dataset"] = adata.obs["dataset"].str.replace("UKIM-V-2", "UKIM-V")

# %%
patients_per_dataset = adata.obs.groupby("dataset").apply(
    lambda x: x["patient"].nunique()
)

# %%
adata_primary_tumor = adata[
    adata.obs["dataset"].str.contains("Lambrechts")
    & (adata.obs["origin_fine"] != "normal_adjacent"),
    :,
]

# %% [markdown]
# ## Infiltration patterns

# %%
ad_immune = sc.AnnData(
    X=(
        adata_primary_tumor.obs.loc[
            lambda x: ~x["cell_type_major"].isin(EXCLUDE_CELL_TYPES)
        ]
        .groupby(
            ["dataset", "patient", "origin_fine", "cell_type_major"], observed=True
        )
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type_major",
            index=["dataset", "patient", "origin_fine"],
            fill_value=0,
        )
    )
)
ad_immune.obs = ad_immune.obs.reset_index().set_index("patient")

# %%
sc.pp.normalize_total(ad_immune, target_sum=1)

# %%
sc.pl.matrixplot(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby=["origin_fine", "patient"],
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
sc.pp.regress_out(ad_immune, "dataset")

# %%
# # for Leader_Merad, one patient has been measured in multiple datasets
# ad_immune = sh.util.aggregate_duplicate_obs(ad_immune)
# assert ad_immune.obs_names.is_unique

# %%
ad_immune.obs["patient"] = ad_immune.obs.index.astype("category")
ad_immune.obs.index.name = "index"

# %%
sc.tl.dendrogram(ad_immune, groupby="patient", use_rep="X", optimal_ordering=True)

# %%
sc.pl.matrixplot(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby="patient",
    dendrogram=True,
    swap_axes=True,
    cmap="bwr",
    vmin=-0.5,
    vmax=0.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %%
sc.pp.neighbors(ad_immune, use_rep="X", n_neighbors=15, metric="correlation")

# %%
sc.tl.leiden(ad_immune, resolution=1)

# %%
sc.pl.heatmap(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby="leiden",
    swap_axes=True,
    cmap="bwr",
    vmin=-0.5,
    vmax=0.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %%
alt.Chart(
    ad_immune.obs.groupby(["leiden", "dataset"]).size().reset_index(name="n")
).mark_bar().encode(
    x="dataset",
    y="n",
    color=alt.Color("leiden", scale=sh.colors.altair_scale("leiden")),
)

# %%
ad_immune.obs["immune_type"] = [
    {
        "0": "T",
        "1": "M",
        "2": "desert",
        "3": "B"
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
    vmin=-0.5,
    vmax=0.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %%
ad_immune.obs

# %%
alt.Chart(ad_immune.obs).encode(
    x=alt.X("patient", axis=alt.Axis(labelLimit=500)),
    y="origin_fine",
    color=alt.Color("immune_type", scale=sh.colors.altair_scale("immune_infiltration")),
).mark_rect()

# %% [markdown]
# ## Make figure

# %%
plot_df = (
    adata.obs.loc[
        :,
        [
            "patient",
            "sex",
            "uicc_stage",
            "ever_smoker",
            "age",
            "tumor_stage",
            "study",
            "platform",
        ],
    ]
    .drop_duplicates()
    .set_index("patient")
    .join(ad_immune.obs, how="inner")
)
plot_df["tumor_stage"] = [
    {"early": "early (I/II)", "advanced": "advanced (III/IV)"}.get(x, x)
    for x in plot_df["tumor_stage"].astype(str).str.replace("None", "unknown")
]
plot_df["immune_type"] = plot_df["immune_type"].astype(str)
remap_tumor_type = {
    "PPC": "NOS",
    "LSCC": "LUSC",
    "LCLC": "NOS",
    "NSCLC": "NOS",
}
plot_df["sex"] = plot_df["sex"].cat.add_categories("unknown").fillna("unknown")
plot_df.rename(
    columns={
        "immune_type": "immune_infiltration",
        "stratum": "TMIG",  # tumor micro environment stratum
    },
    inplace=True,
)
# plot_df.sort_values(
#     ["immune_infiltration"],
#     inplace=True,
#     key=lambda x: [{"n/a": -2, "desert": -1}.get(_, _) for _ in x],
# )
del plot_df["patient"]
plot_df.index.name = "patient"
plot_df = plot_df.reset_index()

# %%
# negative control: random assignment of strata
np.random.seed(0)
plot_df["random_stratum"] = np.array(["desert", "M", "T", "mixed"])[
    np.random.randint(0, 4, size=plot_df.shape[0])
]

# %%
plot_df


# %%
def get_row(col, color_scale=None):
    if color_scale is None:
        color_scale = col

    if color_scale is False:
        altair_scale = alt.Undefined
    else:
        altair_scale = (
            sh.colors.altair_scale(color_scale, data=plot_df, data_col=col)
            if color_scale != "tumor_type"
            else sh.colors.altair_scale(color_scale)
        )
    return (
        alt.Chart(
            plot_df.assign(
                ylab=col,
                label=[
                    f"{p}_{o}"
                    for p, o in zip(ad_immune.obs_names, ad_immune.obs["origin_fine"])
                ],
            )
        )
        .mark_rect()
        .encode(
            x=alt.X(
                "label:N",
                axis=alt.Axis(labels=False, ticks=False, title=None),
            ),
            y=alt.Y("ylab", axis=alt.Axis(title=None)),
            color=alt.Color(
                col,
                scale=altair_scale,
                legend=alt.Legend(columns=3),
            ),
        )
        .properties(width=800)
    )


p0 = (
    alt.vconcat(
        get_row("sex"),
        get_row("tumor_stage", "tumor_stage_verbose"),
        # get_row("study"),
        # get_row("infiltration_state"),
        # get_row("immune_infiltration"),
        get_row("immune_infiltration"),
        get_row("patient", False),
        get_row("origin_fine", False)
    ).resolve_scale(color="independent")
    # .resolve_legend("shared")
    .configure_concat(spacing=0)
)
p0

# %%
heatmap_df = (
    pd.DataFrame(
        ad_immune.X,
        columns=ad_immune.var_names,
        index=[
            f"{p}_{o}"
            for p, o in zip(ad_immune.obs_names, ad_immune.obs["origin_fine"])
        ],
    )
    .fillna(0)
    .reset_index()
    .rename(columns={"index": "label"})
    .melt(id_vars="label")
)
p2 = (
    alt.Chart(heatmap_df)
    .mark_rect()
    .encode(
        x=alt.X(
            "label:N",
            axis=alt.Axis(ticks=True, labels=True, labelLimit=2000),
        ),
        y=alt.Y(
            "cell_type_major",
            sort=[
                "Tumor cells",
                "B cell",
                "Plasma cell",
                "Mast cell",
                "Macrophage",
                "Macrophage FABP4+",
                "Monocyte",
                "DC mature",
                "cDC1",
                "cDC2",
                "pDC",
                "T cell CD4",
                "T cell CD8",
                "T cell regulatory",
                "NK cell",
            ],
            axis=alt.Axis(title=None),
        ),
        color=alt.Color(
            "value",
            scale=sh.colors.altair_scale_mpl("bwr", reverse=False, domain=[-0.5, 0.5]),
        ),
    )
    .properties(width=800, height=200)
)

# %%
(p0 & p2).resolve_scale(x="shared")

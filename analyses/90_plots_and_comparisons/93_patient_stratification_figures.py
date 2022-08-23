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
import statsmodels.api as sm
from tqdm.contrib.concurrent import process_map
import itertools
from operator import or_
from functools import reduce
import scanpy_helpers as sh

# %%
patient_stratification_path = nxfvars.get(
    "patient_stratification_path",
    "../../data/30_downstream_analyses/stratify_patients/stratification/artifacts/patient_stratification.csv",
)
ad_immune_path = nxfvars.get(
    "ad_immune_path",
    "../../data/30_downstream_analyses/stratify_patients/stratification/artifacts/adata_immune.h5ad",
)
ad_tumor_subtypes_path = nxfvars.get(
    "ad_tumor_subtypes_path",
    "../../data/30_downstream_analyses/stratify_patients/stratification/artifacts/adata_tumor_subtypes.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
plot_df = pd.read_csv(patient_stratification_path)
ad_immune = sc.read_h5ad(ad_immune_path)
ad_tumor_subtypes = sc.read_h5ad(ad_tumor_subtypes_path)

# %%
plot_df["immune_infiltration"].value_counts()

# %%
plot_df = plot_df.assign(tumor_type_transcriptomic = lambda x: x["tumor_type_inferred"].map({"LUAD": "LUAD", "LUAD NE": "LUAD", "NSCLC mixed": "NSCLC NOS", "LUAD EMT": "LUAD", "LUSC": "LUSC"}))


# %%
def get_row(col, color_scale=None, title=None, legend_title=None):
    title = col if title is None else title
    legend_title = title if legend_title is None else legend_title
    if color_scale is None:
        color_scale = col
    return (
        alt.Chart(plot_df.assign(ylab=title))
        .mark_rect()
        .encode(
            x=alt.X(
                "patient:N",
                axis=alt.Axis(labels=False, ticks=False, title=None),
                sort=plot_df["patient"].values,
            ),
            y=alt.Y("ylab", axis=alt.Axis(title=None)),
            color=alt.Color(
                col,
                scale=sh.colors.altair_scale(color_scale, data=plot_df, data_col=col),
                legend=alt.Legend(columns=3, title=legend_title),
            ),
        )
        .properties(width=800)
    )


p0 = (
    alt.vconcat(
        get_row("tumor_type_annotated", "tumor_type", title="tumor type (histopathological)", legend_title="tumor type") &
        get_row("tumor_type_transcriptomic", "tumor_type", title="tumor type (transcriptomic)", legend_title="tumor type"),
        get_row("sex"),
        get_row("tumor_stage", "tumor_stage_verbose", title="tumor stage"),
        # get_row("study"),
        get_row("platform"),
        # get_row("infiltration_state"),
        # get_row("immune_infiltration"),
        get_row("immune_infiltration", title="immune infiltration"),
    ).resolve_scale(color="independent")
    # .resolve_legend("shared")
    .configure_concat(spacing=0)
)
p0

# %%
tmp_ad = ad_immune[
    plot_df.loc[plot_df["patient"].isin(ad_immune.obs_names), "patient"], :
].copy()
heatmap_df = (
    pd.DataFrame(tmp_ad.X, columns=tmp_ad.var_names, index=tmp_ad.obs_names)
    .reindex(plot_df["patient"])
    .fillna(0)
    .reset_index()
    .rename(columns={"index": "patient"})
    .melt(id_vars="patient")
)
heatmap_df["cell_type_major"] = heatmap_df["cell_type_major"].str.replace("Tumor cells", "Cancer cells")
p2 = (
    alt.Chart(heatmap_df)
    .mark_rect()
    .encode(
        x=alt.X(
            "patient:N",
            sort=plot_df["patient"].values,
            axis=alt.Axis(ticks=False, labels=False),
        ),
        y=alt.Y(
            "cell_type_major",
            sort=[
                "Cancer cells",
                "B cell",
                "Plasma cell",
                "Mast cell",
                "Macrophage",
                "Macrophage alveolar",
                "Monocyte",
                "T cell CD4",
                "T cell CD8",
                "T cell regulatory",
                "NK cell",
                "DC mature",
                "cDC1",
                "cDC2",
                "pDC",
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

# %%
plot_df.query("tumor_stage == 'unknown'")

# %% [markdown]
# ### tumor cluster fractions

# %%
alt.vconcat(
    get_row("tumor_type_annotated", "tumor_type", title="histopahological").properties(width=1600),
    (
        ad_tumor_subtypes.to_df()
        .reindex(plot_df["patient"])
        .reset_index()
        .melt(id_vars="patient", value_name="fraction")
        .pipe(
            lambda x: alt.Chart(x)
            .mark_bar()
            .encode(
                x=alt.X("patient", sort=plot_df["patient"].tolist(), axis=alt.Axis(labelLimit=500)),
                y=alt.Y("fraction", scale=alt.Scale(domain=[0, 1])),
                color=alt.Color(
                    "cell_type_tumor", scale=sh.colors.altair_scale("tumor_type", data=x, data_col="cell_type_tumor"), legend=alt.Legend(title="transcriptomic")
                ),
            )
        )
    ).properties(height=150),
).resolve_scale(x="shared", color="independent")

# %% [markdown]
# ## groups by histological subtype

# %%
tumor_type_total = (
    plot_df.groupby(["tumor_type_annotated"]).size().reset_index(name="total")
)

# %%
tmp_df = (
    plot_df.groupby(["immune_infiltration", "tumor_type_annotated"])
    .size()
    .reset_index(name="n")
    .merge(tumor_type_total)
    .assign(frac_of_total=lambda x: x["n"] / x["total"])
)

# %%
alt.Chart(tmp_df).mark_bar().encode(
    x=alt.X("immune_infiltration", title=None),
    y="frac_of_total",
    color=alt.Color(
        "immune_infiltration", scale=sh.colors.altair_scale("immune_infiltration")
    ),
).facet(column="tumor_type_annotated")

# %%
alt.Chart(tmp_df).mark_bar().encode(
    x=alt.X("tumor_type_annotated", title=None),
    y="n",
    color=alt.Color("tumor_type_annotated", scale=sh.colors.altair_scale("tumor_type")),
).facet(column="immune_infiltration")

# %% [markdown]
# ### Any associations of groups with sex or tumor stage? 

# %%
smf.glm(
    "sex ~ C(immune_infiltration, Treatment('desert')) + dataset",
    data=plot_df.loc[lambda x: x["sex"] != "unknown"],
    family=sm.families.Binomial(),
).fit().summary()

# %%
plot_df.loc[lambda x: x["sex"] != "unknown"]["tumor_stage"].value_counts()

# %%
smf.glm(
    "tumor_stage ~ C(immune_infiltration, Treatment('desert')) + dataset",
    data=plot_df.loc[lambda x: x["tumor_stage"] != "unknown"],
    family=sm.families.Binomial(),
).fit().summary()

# %%

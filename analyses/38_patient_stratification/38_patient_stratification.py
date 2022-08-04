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
threadpool_limits(32)

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
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018"])
].copy()

# %% [markdown]
# # Annotate tumor cells based on gene expression

# %%
adata_tumor_cells = adata_primary_tumor[
    adata_primary_tumor.obs["cell_type_major"] == "Tumor cells",
    :,
]

# %%
ad_tumor_subtypes = sc.AnnData(
    X=adata_tumor_cells.obs.groupby(["patient", "cell_type_tumor"], observed=True)
    .size()
    .reset_index(name="n")
    .assign(
        cell_type_tumor=lambda x: x["cell_type_tumor"]
        .str.replace("Tumor cells ", "")
        .str.replace(" mitotic", "")
    )
    .pivot_table(values="n", columns="cell_type_tumor", index="patient", fill_value=0),
)
ad_tumor_subtypes.obs = ad_tumor_subtypes.obs.join(
    adata_tumor_cells.obs.loc[:, ["patient", "condition"]]
    .drop_duplicates()
    .set_index("patient")
)


# %%
sc.pp.normalize_total(ad_tumor_subtypes, target_sum=1)

# %%
ax = sc.pl.heatmap(
    ad_tumor_subtypes,
    groupby="condition",
    var_names=ad_tumor_subtypes.var_names,
    swap_axes=True,
    figsize=(14, 1.5),
    show=False,
)
ax["heatmap_ax"].get_figure().savefig(
    f"{artifact_dir}/heatmap_predominant_tumor_cell_type.pdf", bbox_inches="tight"
)

# %%
ad_tumor_subtypes.obs["predominant_tumor_subtype"] = ad_tumor_subtypes.var_names[
    np.argmax(ad_tumor_subtypes.X, axis=1)
]

# %%
ad_tumor_subtypes.obs

# %% [markdown]
# ## Infiltration patterns

# %%
ad_immune = sc.AnnData(
    X=(
        adata_primary_tumor.obs.loc[
            lambda x: ~x["cell_type_major"].isin(EXCLUDE_CELL_TYPES)
            & x["patient"].isin(ad_tumor_subtypes.obs.index)
        ]
        .groupby(["dataset", "patient", "cell_type_major"], observed=True)
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type_major",
            index=["dataset", "patient"],
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
    groupby="patient",
    dendrogram=False,
    swap_axes=True,
    cmap="viridis",
    standard_scale="var",
)

# %%
sc.pp.regress_out(ad_immune, "dataset")

# %%
# for Leader_Merad, one patient has been measured in multiple datasets
ad_immune = sh.util.aggregate_duplicate_obs(ad_immune)
assert ad_immune.obs_names.is_unique

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
sc.tl.leiden(ad_immune, resolution=.75)

# %%
sc.pl.heatmap(
    ad_immune,
    var_names=ad_immune.var_names,
    groupby="leiden",
    swap_axes=True,
    cmap="bwr",
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
        "0": "desert",
        "1": "M",
        "2": "T",
        "3": "T",
        "4": "B",
        "5": "B",
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

# %% [markdown]
# ## Make figure

# %%
# Patients with "primary tumor" samples
adata_primary_tumor.obs["patient"].nunique()

# %%
# Minus patients with no "tumor cells" in primary tumor samples (they can't be used to infer tumor type)
adata_tumor_cells.obs["patient"].nunique()

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
    .join(ad_tumor_subtypes.obs, how="left")
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
plot_df["condition"] = [remap_tumor_type.get(x, x) for x in plot_df["condition"]]
plot_df["sex"] = plot_df["sex"].cat.add_categories("unknown").fillna("unknown")
plot_df.rename(
    columns={
        "condition": "tumor_type_annotated",
        "predominant_tumor_subtype": "tumor_type_inferred",
        "immune_type": "immune_infiltration",
        "stratum": "TMIG",  # tumor micro environment stratum
    },
    inplace=True,
)
plot_df.sort_values(
    ["immune_infiltration", "tumor_type_inferred"],
    inplace=True,
    key=lambda x: [{"n/a": -2, "desert": -1}.get(_, _) for _ in x],
)
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

# %% [markdown]
# # Write output

# %%
plot_df.to_csv("{}/patient_stratification.csv".format(artifact_dir))

# %%
ad_immune.write_h5ad(f"{artifact_dir}/adata_immune.h5ad")

# %%
ad_tumor_subtypes.write_h5ad(f"{artifact_dir}/adata_tumor_subtypes.h5ad")

# %%

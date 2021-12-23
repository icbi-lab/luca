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
from hierarchical_bootstrapping.util import gini_index

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
import hierarchical_bootstrapping as hb
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
import muon as mu

# %%
threadpool_limits(32)

# %%
ah = AnnotationHelper()

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
adata_epi = sc.read_h5ad(
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/33_cell_types_epi/artifacts/adata_epithelial.h5ad"
)

# %%
adata_tumor = sc.read_h5ad(
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/33_cell_types_epi/artifacts/adata_tumor.h5ad"
)

# %%
sc.pl.umap(adata, color=["cell_type", "origin"], ncols=1)

# %%
adata.obs["cell_type"].value_counts()

# %%
cell_types = {
    "tumor": set(adata.obs["cell_type_tumor"]) - set(adata.obs["cell_type"]),
    "healthy epithelial": [
        "Alevolar cell type 2",
        "Club",
        "Ciliated",
        "Alevolar cell type 1",
        "Goblet",
    ],
    "immune": [
        "Macrophage FABP4+",
        "T cell CD4",
        "T cell CD8",
        "Macrophage",
        "Monocyte",
        "NK cell",
        "B cell",
        "T cell regulatory",
        "cDC2",
        "Plasma cell",
        "Mast cell",
        # "Granulocytes",
        "DC mature",
        "pDC",
    ],
    "structural": [
        "Endothelial cell",
        "Fibroblast",
        "Fibroblast adventitial",
        "Fibroblast alevolar",
        "Smooth muscle cell",
        "Pericyte",
    ],
}

# %%
adata_primary_tumor = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018_NSCLC", "Maier_Merad_2020_NSCLC"]),
    :,
]

# %%
adata_tumor_cells = adata_primary_tumor[
    adata_primary_tumor.obs["cell_type"] == "Tumor cells",
    :,
]

# %% [markdown]
# ## Tumor subtypes

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
    .pivot_table(values="n", columns="cell_type_tumor", index="patient", fill_value=0)
)
# ad_tumor_subtypes.obs["patient"] = ad_tumor_subtypes.obs_names

# %%
sc.pp.normalize_total(ad_tumor_subtypes, target_sum=1)

# %%
sc.pl.matrixplot(
    ad_tumor_subtypes,
    groupby="patient",
    var_names=ad_tumor_subtypes.var_names,
    swap_axes=True,
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
def get_cell_type_group(ct):
    for group, cts in cell_types.items():
        if ct in cts:
            return group
    return "other"


# %%
major_cell_types_df = (
    (
        adata_primary_tumor.obs.assign(
            cell_type_group=lambda x: x["cell_type_tumor"].apply(get_cell_type_group)
        )
        .groupby(["dataset", "patient", "cell_type_group"], observed=True)
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type_group",
            index=["dataset", "patient"],
            fill_value=0,
        )
        .drop("other", axis="columns")
    )
    .assign(immune_tumor_ratio=lambda x: np.log2((x["immune"] + 1) / (x["tumor"] + 1)))
    .assign(
        structural_tumor_ratio=lambda x: np.log2(
            (x["structural"] + 1) / (x["tumor"] + 1)
        )
    )
)

# %%
major_cell_types_df

# %%
ad_ti_ratio = sc.AnnData(
    obs=major_cell_types_df.drop(
        columns=["immune_tumor_ratio", "structural_tumor_ratio"]
    ),
    X=major_cell_types_df.loc[:, ["immune_tumor_ratio", "structural_tumor_ratio"]],
)

# %%
ad_ti_ratio.obs = ad_ti_ratio.obs.reset_index().set_index("patient")

# %%
sc.pl.matrixplot(
    ad_ti_ratio,
    groupby="patient",
    var_names=ad_ti_ratio.var_names,
    swap_axes=True,
    cmap="bwr",
    vmin=-7,
    vmax=7
    # vcenter=0
)

# %%
sc.pp.regress_out(ad_ti_ratio, "dataset")

# %%
ad_ti_ratio.obs.index.name = "index"
ad_ti_ratio.obs["patient"] = ad_ti_ratio.obs.index

# %%
sc.pl.matrixplot(
    ad_ti_ratio,
    groupby="patient",
    var_names=ad_ti_ratio.var_names,
    swap_axes=True,
    cmap="bwr",
    vmin=-7,
    vmax=7,
    # vcenter=0
)

# %% [markdown]
# ## All cell-types

# %%
ad_cts = sc.AnnData(
    X=(
        adata_primary_tumor.obs.assign(
            cell_type_group=lambda x: x["cell_type_tumor"].apply(get_cell_type_group)
        )
        .loc[lambda x: x["cell_type_group"] != "other", :]
        .groupby(["dataset", "patient", "cell_type"], observed=True)
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type",
            index=["dataset", "patient"],
            fill_value=0,
        )
    )
)
ad_cts.obs = ad_cts.obs.reset_index().set_index("patient")

# %%
sc.pp.normalize_total(ad_cts, target_sum=1)

# %%
sc.pl.matrixplot(
    ad_cts,
    var_names=ad_cts.var_names,
    groupby="patient",
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
sc.pp.regress_out(ad_cts, "dataset")

# %%
ad_cts.obs["patient"] = ad_cts.obs.index
ad_cts.obs.index.name = "index"

# %%
sc.tl.dendrogram(ad_cts, groupby="patient", use_rep="X", optimal_ordering=True)

# %%
sc.pl.matrixplot(
    ad_cts,
    var_names=ad_cts.var_names,
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

# %% [markdown]
# ## immune cell patterns
# (only immune cells) 

# %%
ad_immune = sc.AnnData(
    X=(
        adata_primary_tumor.obs.assign(
            cell_type_group=lambda x: x["cell_type_tumor"].apply(get_cell_type_group)
        )
        .loc[lambda x: x["cell_type_group"] == "immune", :]
        .groupby(["dataset", "patient", "cell_type"], observed=True)
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type",
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
    # vmin=-0.25,
    # vmax=0.25
    # # vmin=0,
    # vmax=1,
    standard_scale="var",
)

# %%
sc.pp.regress_out(ad_immune, "dataset")

# %%
ad_immune.obs["patient"] = ad_immune.obs.index
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
sc.pp.neighbors(ad_immune)

# %%
sc.tl.umap(ad_immune)

# %%
sc.tl.leiden(ad_immune)

# %%
sc.pl.umap(ad_immune, color=["leiden"] + ad_immune.var_names.tolist(), cmap="bwr")

# %% [markdown]
# ## Pathways

# %%
adata_tumor_pb = hb.tl.pseudobulk(adata_tumor, groupby=["dataset", "patient"])

# %%
sc.pp.normalize_total(adata_tumor_pb, target_sum=1e6)
sc.pp.log1p(adata_tumor_pb)

# %%
model = progeny.load_model(
    organism="Human",  # If working with mouse, set to Mouse
    top=1000,  # For sc we recommend ~1k target genes since there are dropouts
)

# %%
progeny.run(
    adata_tumor_pb,  # Data to use
    model,  # PROGENy network
    center=True,  # Center gene expression by mean per cell
    num_perm=1000,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=False,  # Scale values per feature so that values can be compared across cells
    use_raw=False,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # Pathways with less than 5 targets will be ignored
)

# %%
adata_tumor_progeny = progeny.extract(adata_tumor_pb)

# %%
sc.pp.regress_out(adata_tumor_progeny, "dataset")

# %%
adata_tumor_progeny.obs["patient"] = pd.Categorical(adata_tumor_progeny.obs["patient"])

# %%
sc.tl.dendrogram(
    adata_tumor_progeny, groupby="patient", use_rep="X", optimal_ordering=True
)

# %%
sc.pl.matrixplot(
    adata_tumor_progeny,
    var_names=adata_tumor_progeny.var_names,
    groupby="patient",
    cmap="bwr",
    vmin=-2,
    vmax=2,
    swap_axes=True,
    dendrogram=True,
)

# %%
sc.pp.neighbors(adata_tumor_progeny)

# %%
sc.tl.umap(adata_tumor_progeny)

# %%
sc.pl.umap(
    adata_tumor_progeny, color=["dataset"] + adata_tumor_progeny.var_names.tolist()
)

# %% [markdown]
# ### Multi-Omics

# %%
ad_tumor_subtypes.obs

# %%
ad_ti_ratio.obs

# %%
ad_immune.obs

# %%
adata_tumor_progeny.obs_names = adata_tumor_progeny.obs["patient"]

# %%
mdata = mu.MuData(
    {
        # "cancer type": ad_tumor_subtypes,
        "infiltration patterns": ad_ti_ratio,
        "immune": ad_immune,
        # 'pw': adata_tumor_progeny
    }
)

# %%
mu.pp.intersect_obs(mdata)

# %%
for idx, ad in mdata.mod.items():
    # sc.tl.pca(ad)
    sc.pp.neighbors(ad)

# %%
mu.tl.leiden(mdata, resolution=1.5)

# %%
mdata.obs["leiden"]

# %%
ad_tumor_subtypes.obs["leiden"] = mdata.obs["leiden"]

# %%
sc.pl.heatmap(
    ad_tumor_subtypes,
    groupby="leiden",
    var_names=ad_tumor_subtypes.var_names,
    swap_axes=True,
)

# %%
adata_tumor_progeny.obs["leiden"] = mdata.obs["leiden"]

# %%
sc.pl.heatmap(
    adata_tumor_progeny,
    var_names=adata_tumor_progeny.var_names,
    groupby="leiden",
    cmap="bwr",
    vmin=-2,
    vmax=2,
    swap_axes=True,
)

# %%
ad_ti_ratio.obs["leiden"] = mdata.obs["leiden"]

# %%
sc.pp.neighbors(ad_ti_ratio, n_neighbors=25, metric="cosine")

# %%
sc.tl.leiden(ad_ti_ratio, resolution=0.25)

# %%
sc.pl.heatmap(
    ad_ti_ratio,
    groupby="leiden",
    var_names=ad_ti_ratio.var_names,
    swap_axes=True,
    cmap="bwr",
    vmin=-7,
    vmax=7
    # vcenter=0
)

# %%
imm = ad_ti_ratio[:, "immune_tumor_ratio"].X > 0
structural = ad_ti_ratio[:, "structural_tumor_ratio"].X > 0


def get_state(i, s):
    res = []
    res.append("I+" if i else "I-")
    res.append("S+" if s else "S-")
    return "".join(res)


ad_ti_ratio.obs["group"] = [get_state(i, s) for i, s in zip(imm, structural)]

# %%
sc.pl.heatmap(
    ad_ti_ratio,
    groupby="group",
    var_names=ad_ti_ratio.var_names,
    swap_axes=True,
    cmap="bwr",
    vmin=-7,
    vmax=7
    # vcenter=0
)

# %%
ad_ti_ratio.obs["group"].str.contains("I+", regex=False)

# %%
ad_immune.obs_names = ad_immune.obs_names.values.astype(str)

# %%
ad_immune_sub = ad_immune[
    list(
        set(
            ad_ti_ratio[
                ad_ti_ratio.obs["group"].str.contains("I+", regex=False), :
            ].obs_names
        )
        & set(ad_immune.obs_names)
    ),
    :,
]

# %%
sc.pp.neighbors(
    ad_immune_sub,
    metric="cosine",
    n_neighbors=10,
)

# %%
sc.tl.leiden(ad_immune_sub, resolution=0.5)

# %%
# ad_immune.obs["leiden"] = mdata.obs["leiden"]

# %%
sc.pl.heatmap(
    ad_immune_sub,
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

# %%
ad_immune_sub[ad_immune_sub.obs["leiden"] == "4", :].obs

# %%
ad_ti_ratio.obs["immune_type"] = ad_immune_sub.obs["leiden"]
ad_ti_ratio.obs["immune_type"] = [{"0": "T", "1": "M", "2": "B"}.get(x, "") for x in ad_ti_ratio.obs["immune_type"]]

# %%

# %%
ad_immune2 = ad_immune[ad_tumor_subtypes.obs.index, :]
ad_immune2.obs["immune_type"] =  ad_ti_ratio.obs["immune_type"] 
ad_immune2.obs["group"] =  ad_ti_ratio.obs["group"] 
ad_tumor_subtypes.obs["immune_type"] =  ad_ti_ratio.obs["immune_type"] 
ad_tumor_subtypes.obs["group"] =  ad_ti_ratio.obs["group"] 

# %%
sc.pl.heatmap(
    ad_tumor_subtypes,
    groupby=["group", "immune_type"],
    var_names=ad_tumor_subtypes.var_names,
    swap_axes=True,
)
sc.pl.heatmap(
    ad_ti_ratio[ad_tumor_subtypes.obs.index, :],
    groupby=["group", "immune_type"],
    var_names=ad_ti_ratio.var_names,
    swap_axes=True,
    cmap="bwr",
    vmin=-7,
    vmax=7,
    # vcenter=0
)
sc.pl.heatmap(
    ad_immune2,
    var_names=ad_immune.var_names,
    groupby=["group", "immune_type"],
    swap_axes=True,
    cmap="bwr",
    vmin=-0.5,
    vmax=0.5,
    # # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %%

# %%

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
from threadpoolctl import threadpool_limits

# %%
alt.data_transformers.disable_max_rows()

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
main_adata = nxfvars.get(
    "main_adata",
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)
epithelial_adata = nxfvars.get(
    "epithelial_adata",
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/epithelial_cells_annotated.h5ad",
)
tumor_adata = nxfvars.get(
    "tumor_adata",
    "../../data/20_build_atlas/annotate_datasets/33_cell_types_epi/artifacts/adata_tumor.h5ad",
)
threadpool_limits(nxfvars.get("cpus", 16))

# %%
adata = sc.read_h5ad(main_adata)

# %%
adata_epi = sc.read_h5ad(epithelial_adata)

# %%
adata_tumor = sc.read_h5ad(tumor_adata)

# %% [markdown]
# # Stats and metadata

# %%
adata.X.data.size

# %%
adata.shape[0]

# %%
adata.obs["sample"].nunique()

# %%
adata.obs["platform"].nunique()

# %%
adata.obs["cell_type"].nunique()

# %%
adata.obs["patient"].nunique()

# %%
adata.obs["dataset"].nunique()

# %%
adata.obs["study"].nunique()

# %%
adata.obs.columns

# %%
adata.obs.loc[lambda x: x["origin"].str.contains("tumor")]["patient"].nunique()

# %%
adata.obs.loc[:, ["patient", "origin"]].drop_duplicates().assign(count=1).pivot_table(
    index="patient", columns="origin", values="count", aggfunc=sum
).to_csv(f"{artifact_dir}/origin_count.csv")

# %%
## Export patient table
patient_metadata = (
    adata.obs.loc[
        :,
        [
            "study",
            "dataset",
            "patient",
            "uicc_stage",
            "tumor_stage",
            "sex",
            "ever_smoker",
            "driver_genes",
            "condition",
            "age",
            "platform",
            "platform_fine",
        ]
        + [x for x in adata.obs.columns if "mutation" in x],
    ]
    .drop_duplicates()
    .sort_values(
        [
            "study",
            "dataset",
            "patient",
        ]
    )
    .reset_index(drop=True)
)

patient_metadata.to_csv(f"{artifact_dir}/patient_table.csv")

# %%
pd.set_option("display.max_rows", 20)

# %% tags=[]
patient_metadata

# %% [markdown]
# # UMAP plots
#
# UMAP overview plots of the "core" atlas (without projected datasets)

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adata,
        color="cell_type_coarse",
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_core_atlas.pdf")

# %% [markdown]
# ## Epithelial cells

# %%
# get rid of excluded cells
adata_epi = adata_epi[
    adata_epi.obs_names[adata_epi.obs_names.isin(adata.obs_names)], :
].copy()

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adata_epi,
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=2,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_core_atlas_epithelial.pdf")

# %%
adata_epi.shape[0]

# %%
adata_epi.obs["cell_type"].value_counts()

# %% [markdown]
# ## Tumor cells

# %%
# get rid of excluded cells
adata_tumor = adata_tumor[
    adata_tumor.obs_names[adata_tumor.obs_names.isin(adata.obs_names)], :
]
adata_tumor = adata_tumor[
    adata_tumor.obs["cell_type"].str.contains("Tumor cells"), :
].copy()
adata_tumor.obs["cell_type"] = (
    adata_tumor.obs["cell_type"]
    .str.replace("LSCC", "LUSC")
    .str.replace("Tumor cells ", "")
)

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adata_tumor,
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=3,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_core_atlas_tumor.pdf")

# %%
adata_tumor.shape[0]

# %%
adata_tumor.obs["cell_type"].value_counts()


# %%
def process_subset(mask):
    print(np.sum(mask))
    adata_sub = adata[mask, :].copy()
    sc.pp.neighbors(adata_sub, use_rep="X_scANVI")
    sc.tl.umap(adata_sub)
    return adata_sub


# %%
adatas = {
    label: process_subset(ct)
    for label, ct in {
        "immune": adata.obs["cell_type_coarse"].isin(
            [
                "T cell",
                "NK cell",
                "Neutrophils",
                "Plasma cell",
                "Mast cell",
                "B cell",
                "pDC",
                "cDC",
                "Macrophage/Monocyte",
            ]
        ),
        "structural": adata.obs["cell_type_coarse"].isin(
            ["Endothelial cell", "Stromal"]
        ),
    }.items()
}

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adatas["immune"],
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=6,
        legend_fontoutline=1,
        frameon=False,
        # add_outline=True,
        size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_core_atlas_immune.pdf")

# %%
adatas["immune"].shape[0]

# %%
adatas["immune"].obs["cell_type"].value_counts()

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adatas["structural"],
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=6,
        legend_fontoutline=1,
        frameon=False,
        # add_outline=True,
        size=4,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_core_atlas_structural.pdf")

# %%
adatas["structural"].shape[0]

# %%
adatas["structural"].obs["cell_type"].value_counts()

# %%

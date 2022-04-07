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
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
threadpool_limits(nxfvars.get("cpus", 16))

# %%
adata_core_tumor = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/33_cell_types_epi/artifacts/adata_tumor.h5ad",
)

# %%
adata_core_epi = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/epithelial_cells_annotated.h5ad",
)

# %%
adata = sc.read_h5ad(main_adata)

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
        # size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_extended_atlas.pdf")

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adata,
        color="cell_type_coarse",
        legend_loc=None,
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        # size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_extended_atlas_no_legend.pdf")

# %%
sh.colors.set_scale_anndata(adata, "cell_type_major")
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adata,
        color="cell_type_major",
        groups=[x for x in adata.obs["cell_type_major"].unique() if x != "other"],
        legend_loc="on data",
        legend_fontsize=6,
        legend_fontoutline=2,
        na_in_legend=False,
        frameon=False,
        # add_outline=True,
        # size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_extended_atlas_cell_type_major.pdf")

# %%
adata.obs["study_extended"] = adata.obs.loc[
    lambda x: x["dataset"].str.contains("Leader") | x["dataset"].str.contains("UKIM-V-2"),
    "study",
].astype(str)

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adata,
        color="study_extended",
        frameon=False,
        # add_outline=True,
        size=1,
        return_fig=True,
        title="",
        alpha=0.5
    )
    fig.savefig(f"{artifact_dir}/umap_extended_atlas_projected_data.pdf")


# %% [markdown]
# # Process subsets

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
# Initialize those with the positions from the core atlas.
adatas["epithelial"] = adata[
    adata.obs["cell_type_coarse"].isin(["Epithelial cell"]), :
].copy()
adatas["tumor"] = adata[adata.obs["cell_type_major"] == "Tumor cells", :].copy()

# %%
sh.colors.set_scale_anndata(adatas["tumor"], "cell_type_tumor")

# %%
for core_adata, tmp_adata in zip(
    [adata_core_tumor, adata_core_epi], [adatas[x] for x in ["tumor", "epithelial"]]
):
    # initalize umap with the original coordinates. Missing values (=query) are initialized with random values.
    init_pos_df = pd.DataFrame(
        core_adata.obsm["X_umap"], index=core_adata.obs_names
    ).reindex(tmp_adata.obs_names)
    for col in init_pos_df.columns:
        na_mask = init_pos_df[col].isnull()
        init_pos_df.loc[na_mask, col] = np.random.uniform(
            np.min(init_pos_df[col]), np.max(init_pos_df[col]), size=np.sum(na_mask)
        )
    sc.pp.neighbors(tmp_adata, use_rep="X_scANVI")
    sc.tl.umap(tmp_adata, init_pos=init_pos_df.values)

# %% [markdown]
# # Visualize subsets

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adatas["structural"],
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_extended_atlas_structural.pdf")

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adatas["immune"],
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=6,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=None,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_extended_atlas_immune.pdf")

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adatas["epithelial"],
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_extended_atlas_epithelial.pdf")

# %%
with plt.rc_context({"figure.figsize": (3, 3), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adatas["tumor"],
        color="cell_type_tumor",
        legend_loc="right margin",
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_extended_atlas_tumor.pdf")

# %%
tumor_markers = {
    "Epithelial": ["EPCAM", "B2M"],
    "Alveolar cell type 1": ["AGER", "CLDN18"],
    "Alveolar cell type 2": ["SFTPC", "SFTPB", "SFTPA1"],
    "Ciliated": ["PIFO", "FOXJ1", "HYDIN", "CFAP299"],
    "Club": ["SCGB3A1", "SCGB3A2"],
    "LUAD": ["CD24", "MUC1", "NAPSA", "NKX2-1", "KRT7", "MSLN"],
    "LUSC": ["KRT5", "KRT6A", "TP63", "NTRK2", "SOX2", "KRT17"],
    "NE": ["CHGA", "SYP", "NCAM1", "TUBA1A"],
    "EMT": ["VIM", "SERPINE1", "CDH1"],
    "non-coding": ["MALAT1", "NEAT1"],
    "mitotic": ["MKI67", "TOP2A"],
    "undifferentiated": ["TACSTD2", "AGR2"],
}

# %%
sc.pl.dotplot(adatas["tumor"], groupby="cell_type_tumor", var_names=tumor_markers)

# %%

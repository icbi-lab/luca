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
from nxfvars import nxfvars
import infercnvpy as cnv
from scanpy_helpers.annotation import AnnotationHelper
from scanpy_helpers import de
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from natsort import natsorted
import dorothea
import progeny
import hierarchical_bootstrapping as hb
from natsort import natsorted
import itertools

# %%
ah = AnnotationHelper()

# %%
sc.set_figure_params(figsize=(4, 4))

# %%
input_adata = nxfvars.get(
    "input_adata",
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/31_cell_types_coarse/by_cell_type/adata_cell_type_coarse_epithelial_cell.umap_leiden.h5ad",
)

artifact_dir = nxfvars.get("artifact_dir", "../../data/zz_epi")

# %%
adata = sc.read_h5ad(input_adata)

# %%
adata.obs["leiden"] = adata.obs["leiden_0.50"]

# %%
sc.pl.umap(adata, color="leiden")

# %% [markdown]
# ### Redo UMAP with PAGA

# %%
sc.tl.paga(adata, groups="leiden")

# %%
sc.pl.paga(adata, color="leiden", threshold=0.2)

# %%
sc.tl.umap(adata, init_pos="paga")

# %%
sc.pl.umap(adata, color="dataset")

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    fig = sc.pl.umap(
        adata,
        color=["cell_type", "EPCAM", "dataset", "condition", "origin"],
        cmap="inferno",
        size=2,
        ncols=3,
        return_fig=True,
    )
    fig.savefig(f"{artifact_dir}/overview_epithelial_cluster.pdf", bbox_inches="tight")

# %% [markdown]
# ### Dorothea/progeny

# %%
regulons = dorothea.load_regulons(
    [
        "A",
        "B",
    ],  # Which levels of confidence to use (A most confident, E least confident)
    organism="Human",  # If working with mouse, set to Mouse
)

# %%
from threadpoolctl import threadpool_limits

# %%
threadpool_limits(50)

# %%
dorothea.run(
    adata,  # Data to use
    regulons,  # Dorothea network
    center=True,  # Center gene expression by mean per cell
    num_perm=100,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # TF with less than 5 targets will be ignored
)

# %%
model = progeny.load_model(
    organism="Human",  # If working with mouse, set to Mouse
    top=1000,  # For sc we recommend ~1k target genes since there are dropouts
)

# %%
progeny.run(
    adata,  # Data to use
    model,  # PROGENy network
    center=True,  # Center gene expression by mean per cell
    num_perm=100,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
)

# %%
adata_progeny = progeny.extract(adata)
adata_dorothea = dorothea.extract(adata)

# %%
sc.pl.umap(
    adata_progeny, color=adata_progeny.var_names, cmap="coolwarm", vmin=-2, vmax=2
)

# %%
sc.pl.matrixplot(
    adata_dorothea,
    var_names=adata_dorothea.var_names,
    groupby="leiden",
    cmap="coolwarm",
    vmax=2,
    vmin=-2,
)

# %% [markdown]
# ## Annotation

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    ah.plot_umap(
        adata,
        filter_cell_type=[
            "Ciliated",
            "Alevolar",
            "Basal",
            "Club",
            "Dividing",
            "Goblet",
            "Ionocyte",
            "Mesothelial",
            "Suprabasal",
        ],
        cmap="inferno",
    )

# %%
ah.plot_dotplot(adata)

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
from toolz.functoolz import reduce
from operator import or_

# %%
cell_type_map = {
    "Alevolar cell type 1": [10],
    "Alevolar cell type 2": [0],
    "Ciliated": [7],
    "Club": [1],
}
cell_type_map.update(
    {
        str(k): [k]
        for k in set(map(int, adata.obs["leiden"].values))
        - reduce(or_, map(set, cell_type_map.values()))
    }
)

# %%
ah.annotate_cell_types(adata, cell_type_map)

# %% [markdown]
# ### Subcluster 15, contains alevolar 2

# %%
adata_15 = adata[adata.obs["cell_type"] == "15", :]

# %%
ah.reprocess_adata_subset_scvi(adata_15, leiden_res=0.5)

# %%
sc.pl.umap(adata_15, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
sc.pl.umap(adata_15, color="dataset")

# %%
ah.plot_umap(adata_15, filter_cell_type=["Alev", "Goblet", "Club"])

# %%
ah.plot_dotplot(adata_15)

# %%
ah.annotate_cell_types(
    adata_15,
    {"Alevolar cell type 2": [0, 1, 2, 4], "Alevolar cell type 1": [5], "Club": [3]},
)

# %%
ah.integrate_back(adata, adata_15)

# %% [markdown]
# ### Subcluster 9, contains alevolar 2

# %%
adata_9 = adata[adata.obs["cell_type"] == "9", :]

# %%
ah.reprocess_adata_subset_scvi(adata_9, leiden_res=0.5)

# %%
sc.pl.umap(adata_9, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
sc.pl.umap(adata_9, color=["origin", "dataset"], wspace=0.5)

# %%
ah.plot_umap(adata_9, filter_cell_type=["Alev", "Goblet", "Club", "Epi"])

# %%
ah.annotate_cell_types(
    adata_9,
    {
        "Alevolar cell type 2": [0, 6, 2, 7, 5],
        "Club": [4, 3, 8],
        "ROS1+ normal epithelial": [1],
    },
)

# %%
ah.integrate_back(adata, adata_9)

# %% [markdown]
# ### Subcluster 11, contains goblet

# %%
adata_11 = adata[adata.obs["cell_type"] == "11", :]

# %%
ah.reprocess_adata_subset_scvi(adata_11, leiden_res=0.5)

# %%
sc.pl.umap(adata_11, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
sc.pl.umap(adata_11, color=["origin", "dataset"], wspace=0.5)

# %%
ah.plot_umap(adata_11, filter_cell_type=["Alev", "Goblet", "Club"])

# %%
ah.plot_dotplot(adata_11)

# %%
ah.annotate_cell_types(
    adata_11,
    {
        "Club": [9],
        "Goblet": [7, 3, 1, 0, 4, 5],
        "11-1": [2, 10, 6],
        "11-2": [8],
    },
)

# %%
ah.integrate_back(adata, adata_11)

# %% [markdown]
# ## Find markers for remaining clusters

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(adata, color="cell_type", size=0.6)

# %%
bdata = hb.tl.bootstrap(
    adata, groupby="cell_type", hierarchy=["dataset", "patient"], n=20, use_raw=True
)

# %%
gini_res["expr"].sort_values()

# %%
gini_res = hb.tl.gini(bdata, groupby="cell_type")

# %%
hb.pl.gini_dotplot(adata, gini_res, groupby="cell_type", max_rank=2)

# %%
hb.pl.gini_matrixplot(bdata, gini_res, groupby="cell_type", cmap="Reds", max_rank=2)

# %% [markdown]
# ## Annotate remaining clusters

# %%
nrows = int(np.ceil(adata.obs["cell_type"].nunique() / 4))
fig, axs = plt.subplots(
    nrows, 4, figsize=(16, nrows * 4), gridspec_kw={"wspace": 0.35, "hspace": 0.35}
)
for c, ax in zip(natsorted(adata.obs["cell_type"].unique()), axs.flatten()):
    sc.pl.umap(adata, color="cell_type", groups=[str(c)], size=1, ax=ax, show=False)

# %%
sc.pl.umap(adata, color=["ALB", "GPM6B"], cmap="inferno")

# %%
adata2 = adata.copy()


# %%
cell_type_map

# %%
adata2.obs["cell_type"] = adata.obs["cell_type"]

cell_type_map = {
    "Hemoglobin+": ["11-2"],
    "tumor cell": [2, 3, 8, 6, 5, 4, 16, 17, 18, 14, "11-1"],
    "Alevolar cell type 2": [12, 19, 13],
    "Neuronal cells": [17],
    "Hepatocytes": [20],
}
for ct in set(adata2.obs["cell_type"]) - set(
    map(str, itertools.chain.from_iterable(cell_type_map.values()))
):
    cell_type_map[ct] = cell_type_map.get(ct, []) + [ct]

ah.annotate_cell_types(
    adata2,
    cell_type_map,
    column="cell_type",
)

# %%
adata = adata2

# %% [markdown]
# ## Tumor cells

# %%
adata_tumor = adata[adata.obs["cell_type"] == "tumor cell", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_tumor, leiden_res=0.5)

# %%
# if not removed, paga plotting fails
del adata_tumor.uns["leiden_colors"]

# %%
sc.tl.paga(adata_tumor, "leiden")

# %%
sc.pl.paga(adata_tumor, color="leiden", threshold=0.2)

# %%
sc.tl.umap(adata_tumor, init_pos="paga")

# %%
sc.pl.umap(adata_tumor, color="leiden")

# %%
sc.pl.umap(
    adata_tumor,
    color=[
        "origin",
        "condition",
        "KRT5",
        "KRT14",
        "KRT6A",
        "TP63",
        "CD24",
        "TACSTD2",
        "MUC20",
        "MUC1",
        "NAPSA",
        "NKX2-1",
        "MUC4",
        "CLDN4",
    ],
    cmap="inferno",
)

# %%
sc.pl.umap(
    adata_tumor,
    color=[
        "KRT5",
        "KRT6A",
        "NTRK2",
        "SOX2",
        "MKI67",
        "TOP2A",
        "COL1A1",
        "VIM",
        "ENO2",
        "NCAM1",
        "leiden",
        "dataset",
        "origin",
        "condition",
    ],
    cmap="inferno",
    wspace=0.35,
)

# %%
ah.annotate_cell_types(
    adata_tumor,
    cell_type_map={
        "Tumor cells LSCC mitotic": [1, 17, 14],
        "Tumor cells LSCC": [0],
        "Tumor cells EMT": [3, 11, 13],
        "Tumor cells LUAD mitotic": [8, 5],
        "Tumor cells LUAD": [4, 10, 2, 7, 17, 6, 9, 16, 12],
        "Tumor cells Neuroendocrine": [15],
    },
)

# %%
ah.integrate_back(adata, adata_tumor)

# %%
fig = sc.pl.umap(
    adata,
    color="cell_type",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/marker_umaps.pdf", bbox_inches="tight")

# %% [markdown]
# ### Get markers

# %%
clusters = adata.obs["cell_type"].unique()


def make_means():
    means = np.vstack(
        [np.mean(adata.X[bootstrap(adata, c), :], axis=0).A1 for c in clusters]
    )
    return means


# %%
mean_distribution = np.dstack([make_means() for _ in trange(20)])

# %%
mean_distribution.shape

# %%
mean_distribution.shape

# %%
median_expr = np.median(mean_distribution, axis=2)

# %%
from scipy.stats import rankdata

# %%
gene_ranks = rankdata(-median_expr, axis=0, method="min")

# %%
gini_dist = np.vstack(
    [
        np.apply_along_axis(gini, 0, mean_distribution[:, :, i])
        for i in trange(mean_distribution.shape[2])
    ]
)

# %%
gini_dist.shape

# %%
median_gini = np.median(gini_dist, axis=0)

# %%
gini_df = (
    pd.melt(
        pd.DataFrame(gene_ranks, index=clusters, columns=adata.var_names).reset_index(),
        id_vars=["index"],
        var_name="gene_id",
        value_name="rank",
    )
    .rename(columns={"index": "leiden"})
    .set_index("gene_id")
    .join(pd.DataFrame(index=adata.var_names).assign(gini=median_gini))
    .reset_index()
    .rename(columns={"index": "gene_id"})
    .set_index(["gene_id", "leiden"])
    .join(
        pd.melt(
            pd.DataFrame(
                median_expr, index=clusters, columns=adata.var_names
            ).reset_index(),
            id_vars=["index"],
            var_name="gene_id",
            value_name="expr",
        )
        .rename(columns={"index": "leiden"})
        .set_index(["gene_id", "leiden"])
    )
)

# %%
gini_df.reset_index().sort_values(
    ["leiden", "rank", "gini"], ascending=[True, True, False]
).loc[:, ["leiden", "gene_id", "rank", "gini", "expr"]].to_csv(
    f"{artifact_dir}/markers_all_final.csv"
)

# %%
gini_df.reset_index().sort_values(
    ["leiden", "rank", "gini"], ascending=[True, True, False]
).loc[
    lambda x: (x["expr"] >= 0.5) & (x["rank"] <= 3) & (x["gini"] > 0.6),
    ["leiden", "gene_id", "rank", "gini", "expr"],
].to_csv(
    f"{artifact_dir}/markers_filtered_final.csv"
)

# %%
marker_genes = (
    gini_df.loc[(gini_df["rank"] == 1) & (gini_df["expr"] >= 0.5), :]
    .groupby("leiden")
    .apply(
        lambda df: [
            x[0] for x in df.sort_values("gini", ascending=False).index[:10].values
        ]
    )
    .to_dict()
)

# %%
fig = sc.pl.dotplot(adata, var_names=marker_genes, groupby="cell_type", return_fig=True)
fig.savefig(f"{artifact_dir}/marker_dotplot_final.pdf", bbox_inches="tight")

# %% [markdown]
# ### Write output

# %%
adata.write_h5ad(f"{artifact_dir}/adata_epithelial_cells.h5ad")
adata_tumor.write_h5ad(f"{artifact_dir}/adata_tumor_cells.h5ad")

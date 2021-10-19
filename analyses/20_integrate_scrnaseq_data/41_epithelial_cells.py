# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

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

# %%
ah = AnnotationHelper()

# %%
sc.set_figure_params(figsize=(4, 4))

# %%
input_adata = nxfvars.get(
    "input_adata",
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/split_adata_cell_type/adata_cell_type_coarse_epithelial_cell.umap_leiden.h5ad",
)

artifact_dir = nxfvars.get("artifact_dir", "../../data/zz_epi")

# %%
adata = sc.read_h5ad(input_adata)

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

# %%
adata.obs["leiden"] = adata.obs["leiden_0.50"]

# %%
nrows = int(np.ceil(adata.obs["leiden"].nunique() / 4))
fig, axs = plt.subplots(
    nrows, 4, figsize=(16, nrows * 4), gridspec_kw={"wspace": 0.35, "hspace": 0.35}
)
for c, ax in zip(sorted(adata.obs["leiden"].unique().astype(int)), axs.flatten()):
    sc.pl.umap(adata, color="leiden", groups=[str(c)], size=1, ax=ax, show=False)

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
cell_type_map = {
    "Alevolar cell type 1": [10],
    "Alevolar cell type 2": [0, 12, 17],
    "Ciliated": [7],
    "Tumor cells (distant metastases)": [15, 13],
}
cell_type_map.update(
    {str(k): [k] for k in [1, 14, 16, 19, 4, 6, 3, 5, 2, 11, 8, 18, 9]}
)

# %%
ah.annotate_cell_types(adata, cell_type_map)

# %% [markdown]
# ### Subclustering of cluster 9 (contains goblet cells) 

# %%
adata_9 = adata[adata.obs["leiden"] == "9", :]

# %%
ah.reprocess_adata_subset_scvi(adata_9, leiden_res=0.2)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata_9, color=["leiden", "origin", "condition"])

# %%
ah.plot_dotplot(adata_9)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    ah.plot_umap(
        adata_9, cmap="inferno", filter_cell_type=["Club", "Goblet", "Suprabasal"]
    )

# %%
ah.annotate_cell_types(adata_9, {"Club/Goblet": [0, 1], "9-1": [2, 3], "9-2": [4]})

# %%
ah.integrate_back(adata, adata_9)

# %% [markdown]
# ### Subcluster cluster 8, contains more alevolar cells

# %%
adata_8 = adata[adata.obs["leiden"] == "8", :]

# %%
ah.reprocess_adata_subset_scvi(adata_8, leiden_res=0.2)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata_8, color=["leiden", "origin", "condition"])

# %%
ah.plot_dotplot(adata_8)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    ah.plot_umap(
        adata_8, cmap="inferno", filter_cell_type=["Club", "Goblet", "Alevolar"]
    )

# %%
ah.annotate_cell_types(adata_8, {"Alevolar cell type 2": [0, 4, 5], "8-1": [1, 2, 3]})

# %%
ah.integrate_back(adata, adata_8)

# %% [markdown]
# ### Subcluster cluster 11, contains more alevolar cells

# %%
adata_11 = adata[adata.obs["leiden"] == "11", :]

# %%
ah.reprocess_adata_subset_scvi(adata_11, leiden_res=0.2)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata_11, color=["leiden", "origin", "condition"])

# %%
ah.plot_dotplot(adata_11)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    ah.plot_umap(
        adata_11, cmap="inferno", filter_cell_type=["Club", "Goblet", "Alevolar"]
    )

# %%
ah.annotate_cell_types(
    adata_11, {"Alevolar cell type 2": [2], "11-1": [0], "11-2": [1]}
)

# %%
ah.integrate_back(adata, adata_11)

# %%
with plt.rc_context({"figure.figsize": (7, 7)}):
    fig = sc.pl.umap(adata, color="cell_type", size=1, return_fig=True)
    fig.savefig(f"{artifact_dir}/umap_cell_type.pdf", bbox_inches="tight")

# %% [markdown]
# ## Hierarchical bootstrapping

# %%
import numpy as np
from collections import Counter
import itertools
from numba import jit, njit
from tqdm import trange


# %%
def bootstrap(adata, c):
    idx = []
    leiden_mask = adata.obs["cell_type"] == c
    datasets = np.random.choice(
        adata.obs["dataset"][leiden_mask].unique(), size=np.sum(leiden_mask)
    )
    for d in np.unique(datasets):
        dataset_mask = (adata.obs["dataset"] == d) & leiden_mask
        patients = np.random.choice(
            adata.obs["patient"][dataset_mask].unique(), size=np.sum(dataset_mask)
        )
        for p, p_count in zip(*np.unique(patients, return_counts=True)):
            patient_mask = (adata.obs["patient"] == p) & dataset_mask
            patient_idx = np.where(patient_mask)[0]
            idx.extend(np.random.choice(patient_idx, size=p_count))
    return idx


# %%
sc.pp.normalize_total(adata)

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
@njit
def gini(array):
    """
    Calculate the Gini coefficient of a numpy array.
    Based on: https://github.com/oliviaguest/gini
    Args:
        array (array-like): input array
    Returns:
        float: gini-index of ``array``
    >>> a = np.zeros((10000))
    >>> a[0] = 1.0
    >>> '%.3f' % gini(a)
    '1.000'
    >>> a = np.ones(100)
    >>> '%.3f' % gini(a)
    '0.000'
    >>> a = np.random.uniform(-1,0,1000000)
    >>> '%.2f' % gini(a)
    '0.33'
    """
    # based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    array += 1e-12  # values cannot be 0
    array = np.sort(array)  # values must be sorted
    index = np.arange(1, array.shape[0] + 1)  # index per array element
    n = array.shape[0]  # number of array elements
    return (np.sum((2 * index - n - 1) * array)) / (
        n * np.sum(array)
    )  # Gini coefficient


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
    f"{artifact_dir}/markers_all.csv"
)

# %%
gini_df.reset_index().sort_values(
    ["leiden", "rank", "gini"], ascending=[True, True, False]
).loc[
    lambda x: (x["expr"] >= 0.5) & (x["rank"] <= 3) & (x["gini"] > 0.6),
    ["leiden", "gene_id", "rank", "gini", "expr"],
].to_csv(
    f"{artifact_dir}/markers_filtered.csv"
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
marker_genes

# %%
fig = sc.pl.dotplot(adata, var_names=marker_genes, groupby="cell_type", return_fig=True)
fig.savefig(f"{artifact_dir}/marker_dotplot.pdf", bbox_inches="tight")

# %%
fig = sc.pl.umap(
    adata,
    color=list(itertools.chain(*marker_genes.values())),
    cmap="inferno",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/marker_umaps.pdf", bbox_inches="tight")

# %%
nrows = int(np.ceil(adata.obs["cell_type"].nunique() / 4))
fig, axs = plt.subplots(
    nrows, 4, figsize=(16, nrows * 4), gridspec_kw={"wspace": 0.5, "hspace": 0.5}
)
for c, ax in zip(natsorted(adata.obs["cell_type"].unique()), axs.flatten()):
    sc.pl.umap(adata, color="cell_type", groups=[str(c)], size=1, ax=ax, show=False)

fig.savefig(f"{artifact_dir}/clusters.pdf", bbox_inches="tight")

# %%
adata.write_h5ad(f"{artifact_dir}/adata_epithelial_cells.h5ad")

# %%

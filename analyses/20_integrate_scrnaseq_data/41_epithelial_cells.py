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

# %%
ah = AnnotationHelper()

# %%
sc.set_figure_params(figsize=(4, 4))

# %%
adata = sc.read_h5ad(
    nxfvars.get(
        "input_adata",
        "../../data/20_integrate_scrnaseq_data/annotate_datasets/split_adata_cell_type/adata_cell_type_coarse_epithelial_cell.umap_leiden.h5ad",
    )
)

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.umap(
        adata,
        color=["cell_type", "EPCAM", "dataset", "condition", "origin"],
        cmap="inferno",
        size=2,
        ncols=3,
    )

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
cell_type_map.update({str(k): [k] for k in [1, 14, 16, 19, 4, 6, 3, 5, 2, 11, 8, 18, 9]})

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
ah.annotate_cell_types(adata_11, {"Alevolar cell type 2": [2], "11-1": [0], "11-2": [1]})

# %%
ah.integrate_back(adata, adata_11)

# %%
with plt.rc_context({"figure.figsize": (7, 7)}):
    sc.pl.umap(adata, color="cell_type", legend_loc="on data", legend_fontoutline=1)

# %% [markdown]
# ### EdgeR paired

# %%

# %%
de_res = {}
for c in adata.obs["leiden"].unique():
    try:
        de_res[c] = pd.read_csv(
            f"/home/sturm/projects/2020/pircher-scrnaseq-lung/data/20_integrate_scrnaseq_data/annotate_datasets/de_epi/edger/adata_cell_type_coarse_epithelial_cell_{c}_for_de_de_res_leiden_0.50_edger.tsv",
            #             f"/home/sturm/projects/2020/pircher-scrnaseq-lung/data/20_integrate_scrnaseq_data/annotate_datasets/de_epi/edger_n_cells/adata_cell_type_coarse_epithelial_cell_{c}_for_de_de_res_leiden_0.50_edger.tsv",
            sep="\t",
        ).assign(leiden=c)
    except IOError:
        pass

# %%
de_res_all = pd.concat(de_res.values())

# %%
de_res_all["score"] = (
    -np.log10(de_res_all["FDR"])
    * (de_res_all["logCPM"] > 5)
    * (de_res_all["logFC"] < 1)
)

# %%
de.de_res_to_anndata(
    adata,
    de_res_all,
    groupby="leiden",
    gene_id_col="gene_id",
    pval_col="PValue",
    pval_adj_col="FDR",
    lfc_col="logFC",
    score_col="score",
)

# %%
de_res_all.sort_values("score", ascending=False)

# %%
sc.pl.rank_genes_groups_dotplot(adata, dendrogram=False)

# %% [markdown]
# ### Mast

# %%
de_res = {}
for c in adata.obs["leiden"].unique():
    try:
        de_res[c] = (
            pd.read_csv(
                f"/home/sturm/projects/2020/pircher-scrnaseq-lung/data/20_integrate_scrnaseq_data/annotate_datasets/de_epi/mast/adata_cell_type_coarse_epithelial_cell_{c}_for_de_de_res_leiden_0.50_mast.tsv",
                sep="\t",
            )
            .sort_values("Pr(>Chisq)")
            .assign(leiden=c)
        )
    except IOError:
        pass

# %%
de_res_all = pd.concat(de_res.values())

# %%
de_res_all["score"] = (de_res_all["Pr(>Chisq)"] == 0) * -de_res_all["coef"]

# %%
de.de_res_to_anndata(
    adata,
    de_res_all,
    groupby="leiden",
    gene_id_col="primerid",
    pval_col="Pr(>Chisq)",
    pval_adj_col="Pr(>Chisq)",
    lfc_col="coef",
    score_col="score",
)

# %%
sc.pl.rank_genes_groups_dotplot(adata, dendrogram=False)

# %%
ah.plot_umap(
    adata,
    filter_cell_type=[
        "Alevolar",
        "Basal",
        "Club",
        "Dividing",
        "Goblet",
        "Ionocyte",
        "Mesothelial",
        "Suprabasal",
    ],
    size=1,
    cmap="inferno",
)

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
    leiden_mask = adata.obs["leiden"] == c
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
clusters = adata.obs["leiden"].unique()


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
sc.pl.dotplot(adata, var_names=marker_genes, groupby="leiden")

# %%
sc.pl.umap(adata, color=list(itertools.chain(*marker_genes.values())), cmap="inferno")

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data")

# %%

# %%

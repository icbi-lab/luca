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
sc.set_figure_params(figsize=(5, 5))

# %%
input_adata_epi = nxfvars.get(
    "input_adata_epi",
    "../../data/zz_epi/adata_epithelial_cells.h5ad",
)
input_adata_all = nxfvars.get(
    "input_adata_all",
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/annotate_cell_types_coarse/artifacts/adata_cell_type_coarse.h5ad",
)

artifact_dir = nxfvars.get("artifact_dir", "../../data/zz_epi_deconv")

# %%
adata = sc.read_h5ad(input_adata_all)

# %%
adata_epi = sc.read_h5ad(input_adata_epi)

# %%
ah = AnnotationHelper()

# %%
ah.integrate_back(adata, adata_epi)

# %%
sc.pl.umap(adata, color="cell_type")

# %%
sc.pl.umap(adata_epi, color="cell_type")

# %%

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
    f"{artifact_dir}/markers_all_preliminary.csv"
)

# %%
gini_df.reset_index().sort_values(
    ["leiden", "rank", "gini"], ascending=[True, True, False]
).loc[
    lambda x: (x["expr"] >= 0.5) & (x["rank"] <= 3) & (x["gini"] > 0.6),
    ["leiden", "gene_id", "rank", "gini", "expr"],
].to_csv(
    f"{artifact_dir}/markers_filtered_preliminary.csv"
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
fig.savefig(f"{artifact_dir}/marker_dotplot_preliminary.pdf", bbox_inches="tight")

# %%

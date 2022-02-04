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
import sys

sys.path.insert(0, "/home/sturm/projects/2021/sc-hierarchical-bootstrapping/")

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
from numba import njit
import itertools

# %%
import hierarchical_bootstrapping as hb

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
input_adata_epi = nxfvars.get(
    "input_adata_epi",
    "../../data/zz_epi/adata_epithelial_cells.h5ad",
)
input_adata_all = nxfvars.get(
    "input_adata_all",
    "../../data/0_backup/2021-11-16/20_integrate_scrnaseq_data/annotate_datasets/annotate_cell_types_coarse/artifacts/adata_cell_type_coarse.h5ad",
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
adata.obs.groupby("dataset").apply(lambda x: x.patient.nunique())

# %%
set1 = {
    "Adams_Kaminski_2020_COPD": 44,
    "Lambrechts_2018_LUAD_6149v2": 3,
    "UKIM-V": 3,
    "Kim_Lee_2020_LUAD": 3,
}
set2 = {
    "Reyfman_Misharin_2018_pulmonary-fibrosis": 3,
    "Habermann_Kropski_2020_pulmonary-fibrosis": 3,
    "Mayr_Schiller_2020_pulmonary-fibrosis": 3,
    "He_Fan_2021_LUAD": 44,
}

# %%
gene = "CXCL13"

# %%
datasets_comp = {**set1, **set2}

# %%
patients_by_dataset = {
    dataset: adata.obs.loc[adata.obs["dataset"] == dataset, "patient"]
    .unique()
    .tolist()[: datasets_comp[dataset]]
    for dataset in datasets_comp.keys()
}
selected_patients = list(itertools.chain.from_iterable(patients_by_dataset.values()))

# %%
adata_comp = adata[
    adata.obs["patient"].isin(selected_patients) & (adata.obs["cell_type"] == "T cell"),
    :,
]
adata_comp.obs["group"] = [
    "set1" if d in set1 else "set2" for d in adata_comp.obs["dataset"]
]

# %%
adata_comp.obs.groupby(["group", "dataset"]).apply(lambda x: x.patient.nunique())

# %%
adata_comp.shape

# %%
sc.pl.umap(adata_comp, color=["group", "dataset"])

# %%
sc.pp.normalize_total(adata_comp)
sc.pp.log1p(adata_comp)

# %%
adata_comp.X[:, adata.var_names == gene].todense().A1

# %%
df = (
    pd.DataFrame()
    .assign(
        gene=adata_comp.X[:, adata.var_names == gene].todense().A1,
        group=adata_comp.obs["group"].values,
        patient=adata_comp.obs["patient"].values,
        dataset=adata_comp.obs["dataset"].values,
    )
    .rename(columns={"gene": gene})
)

# %%
import scipy.stats

# %% [markdown]
# ### default method by cell

# %%
scipy.stats.mannwhitneyu(
    df.loc[df["group"] == "set1", gene], df.loc[df["group"] == "set2", gene]
)

# %%
np.mean(df.loc[df["group"] == "set1", gene]), np.mean(
    df.loc[df["group"] == "set2", gene]
)

# %%
df.boxplot(column=gene, by="group")

# %% [markdown]
# ### by patient

# %%
df

# %%
df_agg_patient = (
    df.groupby(["patient", "group"], observed=True).agg(np.mean).reset_index()
)

# %%
scipy.stats.mannwhitneyu(
    df_agg_patient.loc[df_agg_patient["group"] == "set1", gene],
    df_agg_patient.loc[df_agg_patient["group"] == "set2", gene],
)

# %%
df_agg_patient.groupby("group").apply(lambda x: np.mean(x[gene]))

# %%
df_agg_patient.boxplot(column=gene, by="group")

# %% [markdown]
# ### Pseudobulk with balanced patients

# %%
import itertools

# %%
patients_by_dataset = {
    dataset: df.loc[df["dataset"] == dataset, "patient"].unique().tolist()[:3]
    for dataset in df["dataset"].unique()
}
selected_patients = list(itertools.chain.from_iterable(patients_by_dataset.values()))

# %%
df_balanced_patients = df.loc[df["patient"].isin(selected_patients), :]
df_balanced_patients_agg = (
    df_balanced_patients.groupby(["patient", "group"], observed=True)
    .agg(np.mean)
    .reset_index()
)

# %%
scipy.stats.mannwhitneyu(
    df_balanced_patients_agg.loc[df_balanced_patients_agg["group"] == "set1", gene],
    df_balanced_patients_agg.loc[df_balanced_patients_agg["group"] == "set2", gene],
)

# %%
df_balanced_patients_agg.groupby("group").apply(lambda x: np.mean(x[gene]))

# %%
df_balanced_patients_agg.boxplot(column=gene, by="group")

# %% [markdown]
# ### Hierarchical bootstrapping

# %%
adata = hb.tl.bootstrap(
    adata_comp, groupby="group", hierarchy=["dataset", "patient"], n=500
)

# %%
df_bootstrap = pd.DataFrame().assign(
    expr=adata.X[:, adata_comp.var_names == gene][:, 0], group=adata.obs["group"].values
)

# %%
df_bootstrap.boxplot(by="group")

# %% [markdown]
# This is not a valid way to compute a p-value from bootstrapping. The p-value will just increase by increasing the number of bootstrap samples. 
#
# Need to do either
#  * permutation test
#  * studentized bootstrap
#  * ...
#  
# See, for instance 
#  * https://stats.stackexchange.com/questions/20701/computing-p-value-using-bootstrap-with-r
#  * https://blogs.sas.com/content/iml/2011/11/02/how-to-compute-p-values-for-a-bootstrap-distribution.html

# %%
scipy.stats.mannwhitneyu(
    df_bootstrap.loc[lambda x: x["group"] == "set1", "expr"],
    df_bootstrap.loc[lambda x: x["group"] == "set2", "expr"],
)

# %%
df_bootstrap.boxplot(by="group")

# %%
scipy.stats.mannwhitneyu(
    df_bootstrap.loc[lambda x: x["group"] == "set1", "expr"],
    df_bootstrap.loc[lambda x: x["group"] == "set2", "expr"],
)

# %%
np.mean(df_bootstrap["set1"]), np.mean(df_bootstrap["set2"])

# %%
adata.shape

# %%
adata.obs

# %%
import numpy as np
from collections import Counter
import itertools
from numba import jit, njit
from tqdm import trange


# %%
# This variant does not appropriately control for false negatives (very small pvalue for housekeeping gene)
def bootstrap(adata, c, col="group"):
    idx = []
    leiden_mask = adata.obs[col] == c
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
def bootstrap2(adata, c, col="group", n_patients=3, n_cells=300):
    idx = []
    leiden_mask = adata.obs[col] == c
    uniq_datasets = adata.obs["dataset"][leiden_mask].unique()
    datasets = np.random.choice(uniq_datasets, uniq_datasets.size)
    for d in datasets:
        # print(d)
        dataset_mask = (adata.obs["dataset"] == d) & leiden_mask
        uniq_patients = adata.obs["patient"][dataset_mask].unique()
        patients = np.random.choice(uniq_patients, n_patients)
        for p in patients:
            # print("  " + p)
            patient_mask = (adata.obs["patient"] == p) & dataset_mask
            cells = np.where(patient_mask)[0]
            idx.extend(np.random.choice(cells, n_cells))
    return idx


# %%
# sc.pp.normalize_total(adata)

# %%
adata_comp.obs["group"].unique()


# %%
def make_means():
    means = np.vstack(
        [
            np.mean(adata_comp.X[bootstrap2(adata_comp, c), :], axis=0).A1
            for c in adata_comp.obs["group"].unique()
        ]
    )
    return means


# %%
# %%time
mean_distribution = np.dstack([make_means() for _ in trange(100)])

# %%
adata_comp.shape

# %%
mean_distribution.shape

# %%
mean_distribution.shape

# %% [markdown]
# ---
#
# ## Old gini stuff

# %%
assert False


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

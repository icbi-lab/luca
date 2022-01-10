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
from threadpoolctl import threadpool_limits, threadpool_info

threadpool_limits(32)

# %%
import scanpy as sc
import pandas as pd
from nxfvars import nxfvars

sc.settings.set_figure_params(figsize=(5, 5))
from pathlib import Path
from scanpy_helpers.annotation import AnnotationHelper
import scanpy_helpers as sh
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
import seaborn as sns
import pickle
import warnings
import anndata

warnings.simplefilter(action="ignore", category=FutureWarning)

# TODO remove those filters
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore")

warnings.simplefilter(action="ignore", category=anndata.ImplicitModificationWarning)

# %%
ah = AnnotationHelper()

# %% [markdown]
# # Load input data

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)

# %%
patient_strat = pd.read_csv(
    nxfvars.get(
        "stratification_csv",
        "../../data/30_downstream_analyses/stratify_patients/artifacts/patient_stratification.csv",
    ),
    index_col=0,
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
patient_strat.set_index("patient", inplace=True)
patient_strat["dataset"] = (
    adata.obs.loc[:, ["patient", "dataset"]]
    .drop_duplicates()
    .set_index("patient")["dataset"]
)
patient_strat.reset_index(inplace=True)

# %%
patient_strat

# %% [markdown]
# ## Load squidpy results (cell2cell)
#
# Load the results obtained for each sample and 
# load them into an anndata file (each obs is a sample and each var a pair of receptor/ligand)

# %%
cpdb_res = {}
for f in Path("../../data/30_downstream_analyses/cell2cell/squidpy/").glob("**/*.pkl"):
    sample = f.name.replace("full_atlas_annotated_", "").replace(".pkl", "")
    with open(f, "rb") as fh:
        cpdb_res[sample] = pickle.load(fh)

# %%
dfs_melt = {}
for k in tqdm(cpdb_res):
    dfs_melt[k] = (
        cpdb_res[k]["means"]
        .reset_index()
        .melt(id_vars=["source", "target"], value_name=k)
    )

# %%
var = pd.concat(
    [
        df.loc[lambda x: x[k] != 0, ["source", "target", "cluster_1", "cluster_2"]]
        for k, df in tqdm(dfs_melt.items())
    ]
).drop_duplicates()

# %%
var = var.assign(idx=lambda x: ["_".join(t[1:]) for t in x.itertuples()]).set_index(
    "idx"
)

# %%
for k, df in tqdm(dfs_melt.items()):
    tmp_series = (
        df.loc[lambda x: x[k] != 0, :]
        .assign(idx=lambda x: ["_".join(t[1:-1]) for t in x.itertuples()])
        .set_index("idx")[k]
    )
    var[k] = tmp_series

# %%
ad_cpdb = sc.AnnData(var=var.iloc[:, :4], X=var.iloc[:, 4:].T.fillna(0))

# %%
sample_info = (
    adata.obs.loc[
        :, ["sample", "patient", "tissue", "origin", "condition", "tumor_stage"]
    ]
    .assign(sample_lc=lambda x: x["sample"].str.lower())
    .drop_duplicates()
    .merge(patient_strat, on="patient", how="left")
    .set_index("sample_lc")
)

# %%
ad_cpdb.obs = ad_cpdb.obs.join(sample_info)

# %%
ad_cpdb = ad_cpdb[:, ad_cpdb.var["cluster_1"] != ad_cpdb.var["cluster_2"]].copy()

# %%
ad_cpdb.shape

# %% [markdown]
# # Run differential signature analysis

# %%
cytosig_signature = pd.read_csv(
    "../../tables/cytosig_signature_matrix.tsv", sep="\t", index_col=0
)

# %%
regulons = dorothea.load_regulons(
    [
        "A",
        "B",
    ],  # Which levels of confidence to use (A most confident, E least confident)
    organism="Human",  # If working with mouse, set to Mouse
)

# %%
model = progeny.load_model(
    organism="Human",  # If working with mouse, set to Mouse
    top=1000,  # For sc we recommend ~1k target genes since there are dropouts
)


# %%
@sh.util.suppress_stdout
def run_progeny(adata):
    tmp_adata = adata.copy()
    progeny.run(
        tmp_adata,  # Data to use
        model,  # PROGENy network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # Pathways with less than 5 targets will be ignored
    )
    return progeny.extract(tmp_adata)


# %%
@sh.util.suppress_stdout
def run_dorothea(adata):
    tmp_adata = adata.copy()
    dorothea.run(
        tmp_adata,  # Data to use
        regulons,  # Dorothea network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # TF with less than 5 targets will be ignored
    )
    return dorothea.extract(tmp_adata)


# %%
@sh.util.suppress_stdout
def run_cytosig(adata):
    tmp_adata = adata.copy()
    progeny.run(
        tmp_adata,  # Data to use
        cytosig_signature,  # PROGENy network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # Pathways with less than 5 targets will be ignored
        obsm_key="cytosig",
    )
    return progeny.extract(tmp_adata, obsm_key="cytosig")


# %%
def prepare_cpdb(adata):
    tmp_ad = ad_cpdb[ad_cpdb.obs["sample"].isin(adata.obs["sample"].unique()), :]
    # only keep interactions that a >0 in at least N samples
    tmp_ad = tmp_ad[:, np.sum(tmp_ad.X != 0, axis=0) > 10].copy()
    return tmp_ad


# %% [markdown]
# ## Prepare datasets

# %%
adata.obs["cell_type_structural"] = [
    {
        "Tumor cells": "tumor cells",
        "Alveolar cell type 1": "healthy epithelial",
        "Alveolar cell type 2": "healthy epithelial",
        "Goblet": "healthy epithelial",
        "Club": "healthy epithelial",
        "Ciliated": "healthy epithelial",
        "Fibroblast": "stromal",
        "Fibroblast adventitial": "stromal",
        "Fibroblast alveolar": "stromal",
        "Smooth muscle cell": "stromal",
        "Pericyte": "stromal",
        "Endothelial cell": "endothelial",
    }.get(ct, "other")
    for ct in adata.obs["cell_type"]
]

# %%
adata_tumor_normal = adata[
    adata.obs["origin"].isin(["normal_adjacent", "normal", "tumor_primary"]),
    :,
].copy()

# %%
adata_tumor_normal.obs["origin"] = [
    "normal" if "normal" in x else x for x in adata_tumor_normal.obs["origin"]
]

# %%
adata_tumor_normal.obs["origin"].unique()

# %%
sc.pl.umap(
    adata_tumor_normal,
    color=["cell_type_major"],
    ncols=2,
    wspace=1,
)

# %%
adata_primary_tumor = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"]),
    :,
]

# %%
adata_primary_tumor.obs = (
    adata_primary_tumor.obs.reset_index()
    .merge(
        patient_strat,
        how="left",
        left_on=["patient", "dataset"],
        right_on=["patient", "dataset"],
    )
    .set_index("index")
)

# %%
sc.pl.umap(
    adata_primary_tumor,
    color=["cell_type_major", "cell_type_structural"],
    ncols=2,
    wspace=1,
)

# %%
comparisons = {
    "tumor_normal": {
        "dataset": adata_tumor_normal,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "origin",
        # paired analysis
        "lm_covariate_str": "+ patient",
    },
    "patient_immune_infiltration": {
        "dataset": adata_primary_tumor,
        "cell_type_column": "cell_type_structural",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "immune_infiltration",
        "lm_covariate_str": "+ dataset",
    },
    "luad_lscc": {
        "dataset": adata_primary_tumor[
            adata_primary_tumor.obs["condition"].isin(["LUAD", "LSCC"]), :
        ],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "condition",
        "lm_covariate_str": "+ dataset",
    },
    "early_advanced": {
        "dataset": adata_primary_tumor,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "tumor_stage",
        "lm_covariate_str": "+ dataset",
    },
}

# %% [markdown]
# ---

# %%
TOOLS = {
    "progeny": run_progeny,
    "dorothea": run_dorothea,
    "cytosig": run_cytosig,
}


def prepare_dataset(
    id_, *, dataset, cell_type_column, tools=("progeny", "cytosig"), **kwargs
):
    """Split anndata by cell-type and run the different signature enrichment
    methods"""
    print(f"Preparing dataset for {id_}:")
    dataset = dataset[dataset.obs[cell_type_column] != "other", :]
    print(f"\tSplitting anndata by {cell_type_column}:")
    adata_by_cell_type = sh.util.split_anndata(dataset, cell_type_column)
    all_adatas = {}
    for tool in tools:
        print(f"\tRunning {tool}:")
        all_adatas[tool] = {
            ct: TOOLS[tool](tmp_adata)
            for ct, tmp_adata in tqdm(adata_by_cell_type.items())
        }

    return all_adatas


def compare_signatures(
    id_,
    all_adatas,
    *,
    pseudobulk_group_by,
    column_to_test,
    lm_covariate_str,
    **kwargs,
):
    """Compare signature enrichment using linear model"""
    print(f"Performing comparison for {id_}:")
    all_results = {}
    for tool in all_adatas:
        print(f"\tRunning tests for {tool}:")
        tmp_res = sh.lm.lm_test_all(
            all_adatas[tool],
            groupby=pseudobulk_group_by,
            column_to_test=column_to_test,
            lm_covariate_str=lm_covariate_str,
        )
        tmp_res.to_csv(
            f"{artifact_dir}/differential_signature_{id_}_{tool}.tsv", sep="\t"
        )
        all_results[tool] = tmp_res

    return all_results


def compare_cpdb(
    id_,
    *,
    dataset,
    column_to_test,
    lm_covariate_str,
    **kwargs,
):
    """
    Perform test on cellphonedb samples.
    The data is ALWAYS grouped on sample!
    """
    # filter cpdb anndata object to contain the same samples
    tmp_ad = ad_cpdb[
        ad_cpdb.obs["sample"].isin(dataset.obs["sample"].unique().tolist()), :
    ].copy()
    tmp_res = (
        sh.lm.test_lm(
            tmp_ad,
            f"~ 0 + C({column_to_test}, Sum) {lm_covariate_str}",
            column_to_test,
            progress=True,
            n_jobs=32,
        )
        .set_index("variable")
        .join(ad_cpdb.var, how="left")
        .sort_values("pvalue")
        .dropna(how="any")
        .assign(
            fdr=lambda x: statsmodels.stats.multitest.fdrcorrection(x["pvalue"].values)[
                1
            ]
        )
    )
    tmp_res.to_csv(f"{artifact_dir}/differential_signature_{id_}_cpdb.tsv", sep="\t")
    return tmp_res


# %%
datasets = {id_: prepare_dataset(id_, **config) for id_, config in comparisons.items()}

# %%
results = {
    id_: compare_signatures(id_, datasets[id_], **config)
    for id_, config in comparisons.items()
}

# %%
ad_cpdb.var

# %%
ad_cpdb.obs

# %%
results["tumor_normal"]["cpdb"]

# %%
for id_, config in comparisons.items():
    if "cpdb" in config.get("tools", ["cpdb"]):
        results[id_]["cpdb"] = compare_cpdb(id_, **config)


# %% [markdown]
# ---

# %%
def mk_matrixplot(cell_type, features, adatas, column):
    ad_lookup = {ct: ad for ct, ad in zip(cell_types, adatas)}
    tmp_bulk = hb.tl.pseudobulk(
        ad_lookup[cell_type],
        groupby=["patient", column],
        aggr_fun=np.mean,
    )
    sc.pl.matrixplot(tmp_bulk, var_names=features, groupby=column, title=cell_type)


# %%
for cell_type in res_cytosig["cell_type"].unique():
    mk_matrixplot(
        cell_type,
        res_cytosig.loc[lambda x: x["cell_type"] == cell_type, "variable"],
        adatas_cytosig,
        "immune_infiltration",
    )

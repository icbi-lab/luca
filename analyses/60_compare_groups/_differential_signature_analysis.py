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
import seaborn as sns
import pickle

# %%
ah = AnnotationHelper()

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
# ## Run differential signature analysis

# %%
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
adata.obs["cell_type2"] = [
    {
        "Tumor cells": "tumor cells",
        "Alevolar cell type 1": "healthy epithelial",
        "Alevolar cell type 2": "healthy epithelial",
        "Goblet": "healthy epithelial",
        "Club": "healthy epithelial",
        "Ciliated": "healthy epithelial",
        "Fibroblast": "stromal",
        "Fibroblast adventitial": "stromal",
        "Fibroblast alevolar": "stromal",
        "Smooth muscle cell": "stromal",
        "Pericyte": "stromal",
        "Endothelial cell": "endothelial",
    }.get(ct, "other")
    for ct in adata.obs["cell_type"]
]

# %%
adata_primary_tumor = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018_NSCLC", "Maier_Merad_2020_NSCLC"]),
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
    adata_primary_tumor, color=["cell_type_coarse", "cell_type", "cell_type2"], ncols=1
)

# %%
adatas_by_cell_type = []
cell_types = [x for x in adata_primary_tumor.obs["cell_type2"].unique() if x != "other"]
for cell_type in cell_types:
    adatas_by_cell_type.append(
        adata_primary_tumor[
            adata_primary_tumor.obs["cell_type2"] == cell_type, :
        ].copy()
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
cytosig_signature = pd.read_csv(
    "../../tables/cytosig_signature_matrix.tsv", sep="\t", index_col=0
)

# %%
adatas_dorothea = list(
    tqdm(map(run_dorothea, adatas_by_cell_type), total=len(adatas_by_cell_type))
)

# %%
adatas_progeny = list(
    tqdm(map(run_progeny, adatas_by_cell_type), total=len(adatas_by_cell_type))
)

# %%
adatas_cytosig = list(
    tqdm(map(run_cytosig, adatas_by_cell_type), total=len(adatas_by_cell_type))
)

# %%
all_adatas = {
    "progeny": adatas_progeny,
    "cytosig": adatas_cytosig,
    "dorothea": adatas_dorothea,
}


# %%
def test_lm(pseudobulk, formula, column, progress=False):
    """
    Use a linear model to find differences between groups

    In this case we use sum-to-zero or deviation coding to find
    deviations from the mean of means

    tmp_adata is a pseudobulk anndata object"""
    var_names = pseudobulk.var_names
    
    df = pseudobulk.obs.join(
        pd.DataFrame(pseudobulk.X, index=pseudobulk.obs_names, columns=var_names)
    )

    pseudobulk.obs[column] = pd.Categorical(pseudobulk.obs[column])
    all_groups = pseudobulk.obs[column].unique()
    
    def test_all_params(res, all_groups):
        # only using the categories gets rid of NAN
        keys = [f"C({column}, Sum)[S.{g}]" for g in all_groups.categories.values]
        # print(keys)
        # print(res.params)
        coefs = res.params[keys[:-1]].to_dict()
        pvals = res.pvalues[keys[:-1]].to_dict()
        # test the level that was omitted for redundancy
        coefs[keys[-1]] = -sum(coefs.values())
        pvals[keys[-1]] = float(
            res.f_test(" + ".join([f"{k}" for k in keys[:-1]]) + " = 0").pvalue
        )
        return coefs, pvals

    results = []
    lms = []
    var_iter = tqdm(var_names) if progress else var_names
    for col in var_iter:
        group_results = []
        mod = smf.ols(formula=formula.format(col=col), data=df)
        res = mod.fit()
        coefs, pvals = test_all_params(res, all_groups)
        res_df = (
            pd.DataFrame.from_dict(coefs, orient="index", columns=["coef"])
            .join(pd.DataFrame.from_dict(pvals, orient="index", columns=["pvalue"]))
            .assign(
                variable=col,
                group=lambda x: [
                    re.search("\[S\.(.*)\]", k).groups()[0] for k in x.index
                ],
            )
        )
        results.append(res_df)
        lms.append(res)

    return lms, pd.concat(results)


# %%
res_per_cell_type = {}
column = "immune_infiltration"
for signature, adatas in all_adatas.items():
    tmp_res = []
    for ct, tmp_adata in zip(tqdm(cell_types), adatas):
        tmp_bdata = hb.tl.pseudobulk(
            tmp_adata,
            groupby=["dataset", "patient", column],
            aggr_fun=np.mean,
        )
        if tmp_bdata.obs[column].nunique() < 3:
            continue
        _, res_df = test_lm(
            tmp_bdata, f"Q('{{col}}') ~ 0 + C({column}, Sum) + dataset", column
        )
        res_df = res_df.reset_index(drop=True).assign(cell_type=ct)
        tmp_res.append(res_df)

    res_per_cell_type[signature] = pd.concat(tmp_res).dropna(how="any")

# %%
for tmp_df in res_per_cell_type.values():
    _, tmp_df["fdr"] = statsmodels.stats.multitest.fdrcorrection(
        tmp_df["pvalue"].values
    )

# %%
res_per_cell_type["progeny"].sort_values("pvalue").groupby("group").apply(
    lambda x: x.query("fdr < 0.1")
)

# %%
res_per_cell_type["dorothea"].sort_values("pvalue").groupby("group").apply(
    lambda x: x.query("fdr < 0.1")
)

# %%
res_cytosig = (
    res_per_cell_type["cytosig"]
    .sort_values("pvalue")
    .groupby("group")
    .apply(lambda x: x.query("fdr < 0.1"))
)
res_cytosig


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

# %% [markdown]
# ---
# # CNV

# %%
# TODO: compare to SCEVAN
# TODO: regress out dataset-specific effects, or at least include dataset in linear model.

# %%
ithcna_res_files = list(
    Path("../../data/30_downstream_analyses/infercnv/infercnvpy/").glob("**/ithcna.txt")
)
ithgex_res_files = list(
    Path("../../data/30_downstream_analyses/infercnv/infercnvpy/").glob("**/ithgex.txt")
)
cnvscore_res_files = list(
    Path("../../data/30_downstream_analyses/infercnv/infercnvpy/").glob("**/cnv_score.txt")
)
scevan_res_files = list(
    Path("../../data/30_downstream_analyses/infercnv/scevan/").glob("**/scevan_result.csv")
)
adatas_infercnvpy = list(
    Path("../../data/30_downstream_analyses/infercnv/infercnvpy/").glob("**/*.h5ad")
)

# %%
adatas_cnv = {}
for f in tqdm(adatas_infercnvpy):
    patient = str(f).split("/")[-2].replace("full_atlas_annotated_", "")
    adatas_cnv[patient] = sc.read_h5ad(f)

# %%
cnvsum_res = {}
for dataset, ad in tqdm(adatas_cnv.items()):
    x = ad[ad.obs["cell_type"] == "Tumor cells", :].obsm["X_cnv"].copy()
    x[x < 0] = 0
    cnvsum_res[dataset] = np.mean(np.sum(x, axis=1))

# %%
ithcna_res = {}
for f in ithcna_res_files:
    patient = str(f).split("/")[-2].replace("full_atlas_annotated_", "")
    ithcna = np.loadtxt(f)
    ithcna_res[patient] = float(ithcna)
    
ithgex_res = {}
for f in ithgex_res_files:
    patient = str(f).split("/")[-2].replace("full_atlas_annotated_", "")
    ithgex = np.loadtxt(f)
    ithgex_res[patient] = float(ithgex)
    
cnvscore_res = {}
for f in cnvscore_res_files:
    patient = str(f).split("/")[-2].replace("full_atlas_annotated_", "")
    cnvscore = np.loadtxt(f)
    cnvscore_res[patient] = float(cnvscore)
    
n_subclones = {}
for f in scevan_res_files:
    patient = str(f).split("/")[-2].replace("full_atlas_annotated_", "")
    df = pd.read_csv(f, index_col=0)
    if "subclone" in df.columns:
        n_subclones[patient] = np.max(df["subclone"].dropna())

# %%
patient_strat["patient_lc"] = patient_strat["patient"].str.lower()

# %%
patient_strat.set_index("patient_lc", inplace=True)

# %%
patient_strat["ithcna"] = pd.Series(ithcna_res)
patient_strat["ithgex"] = pd.Series(ithgex_res)
patient_strat["cnvscore"] = pd.Series(cnvscore_res)
patient_strat["n_subclones"] = pd.Series(n_subclones)
patient_strat["cnvsum"] = pd.Series(cnvsum_res)

# %%
pat_t = patient_strat["immune_infiltration"][lambda x: x == "T"].index.values.copy()
pat_none = patient_strat["immune_infiltration"][lambda x: x == "-"].index.values.copy()
np.random.shuffle(pat_t)
np.random.shuffle(pat_none)

# %%
import infercnvpy as cnv

# %%
# TODO use linear model to check for confounding effects: dataset and number of cells! 
for var in ["ithcna", "ithgex", "cnvscore", "n_subclones"]:
# for var in ["cnvsum"]:
    print(var)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(22, 5))
    sns.swarmplot(x="TMIG", y=var, data=patient_strat, ax=ax1, hue="dataset")
    sns.boxplot(x="TMIG", y=var, data=patient_strat, ax=ax1, width=0.2)
    ax1.get_legend().remove()
    
    sns.swarmplot(
        x="infiltration_state", y=var, data=patient_strat, ax=ax2, hue="dataset"
    )
    sns.boxplot(x="infiltration_state", y=var, data=patient_strat, ax=ax2, width=0.2)
    ax2.get_legend().remove()
    
    sns.swarmplot(
        x="tumor_type_inferred", y=var, data=patient_strat, ax=ax3, hue="dataset"
    )
    sns.boxplot(x="tumor_type_inferred", y=var, data=patient_strat, ax=ax3, width=0.2)
    ax3.get_legend().remove()
    
    sns.swarmplot(
        x="immune_infiltration", y=var, data=patient_strat, ax=ax4, hue="dataset"
    )
    sns.boxplot(x="immune_infiltration", y=var, data=patient_strat, ax=ax4, width=0.8)
    ax4.get_legend().remove()
    plt.show()

# %%
scipy.stats.mannwhitneyu(
    patient_strat.loc[lambda x: x["immune_infiltration"] == "-", "cnvscore"],
    patient_strat.loc[lambda x: x["immune_infiltration"] == "B", "cnvscore"],
)

# %% [markdown]
# # Cell2Cell

# %%
cpdb_res = {}
for f in Path("../../data/30_downstream_analyses/cell2cell/squidpy/").glob("*.pkl"):
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
dfs_melt["maynard_bivona_2020_nsclc_lt_s01"]

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
var

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
    adata.obs.loc[:, ["sample", "patient", "tissue", "origin", "condition"]]
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

# %%
# ad_cpdb_nsclc = ad_cpdb[~ad_cpdb.obs["TMIG"].isnull(), :].copy()
ad_cpdb_primary = ad_cpdb[
    adata_primary_tumor.obs["sample"].str.lower().unique(), :
].copy()

# %%
ad_cpdb_primary = ad_cpdb_primary[:, np.sum(ad_cpdb_primary.X != 0, axis=0) > 10].copy()

# %%
ad_cpdb_primary.shape

# %%
ad_cpdb_primary.obs


# %%
def chunk_adatas(ad, chunksize=200):
    for i in range(0, ad.shape[1], chunksize):
        yield ad[:, i:i+chunksize].copy()


# %%
def do_test(adata):
    lms, df = test_lm(
        adata,
        "Q('{col}') ~ 0 + C(TMIG, Sum) + dataset",
        "TMIG",
        progress=True,
    )
    return df


# %%
res_df = pd.concat(process_map(do_test, list(chunk_adatas(ad_cpdb_primary)), max_workers=32))

# %%
_, res_df["fdr"] = statsmodels.stats.multitest.fdrcorrection(res_df["pvalue"].values)

# %%
pd.set_option("display.max_rows", None)

# %%
res_df.query("fdr < 0.1").set_index("variable").join(ad_cpdb_primary.var).sort_values(
    ["group", "pvalue"]
)

# %%
res_df.query("fdr < 0.1").set_index("variable").join(ad_cpdb_primary.var).sort_values(
    ["group", "cluster_1", "cluster_2"]
)

# %%

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
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
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
    sample = f.name.replace("full_atlas_merged_", "").replace(".pkl", "")
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
        :,
        [
            "sample",
            "patient",
            "tissue",
            "origin",
            "condition",
            "tumor_stage",
            "dataset",
        ],
    ]
    .assign(sample_lc=lambda x: x["sample"].str.lower())
    .drop_duplicates()
    .merge(patient_strat, how="left")
    .set_index("sample_lc")
)

# %%
ad_cpdb.obs = ad_cpdb.obs.join(sample_info)

# %%
ad_cpdb = ad_cpdb[
    :,
    (ad_cpdb.var["cluster_1"] != ad_cpdb.var["cluster_2"])
    & (ad_cpdb.var["cluster_1"] != "other")
    & (ad_cpdb.var["cluster_2"] != "other"),
].copy()

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
    )
    .set_index("index")
)

# %%
adata_primary_tumor.obs["infiltration_status"] = [
    {"-": "immune_low", "T": "immune_high", "B": "immune_high", "M": "immune_high"}.get(
        i, np.nan
    )
    for i in adata_primary_tumor.obs["immune_infiltration"]
]
adata_primary_tumor.obs["infiltration_type"] = [
    np.nan if i == "-" else i for i in adata_primary_tumor.obs["immune_infiltration"]
]
adata_primary_tumor.obs["ct_all"] = [
    "all" if ct != "other" else ct for ct in adata_primary_tumor.obs["cell_type_major"]
]

# %%
# TODO: split and aggregate datasets appropriately for cpdb.
comparisons = {
    "tumor_normal": {
        "dataset": adata_tumor_normal,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "origin",
        # paired analysis
        "lm_covariate_str": "+ patient",
        "contrasts": "Treatment('normal')",
        "tools": ["dorothea", "progeny", "cytosig", "cpdb"],
    },
    "infiltration_status": {
        "dataset": adata_primary_tumor[
            adata_primary_tumor.obs["infiltration_status"].isin(
                ["immune_low", "immune_high"]
            ),
            :,
        ],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "infiltration_status",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Treatment('immune_low')",
    },
    "infiltration_type": {
        "dataset": adata_primary_tumor[
            adata_primary_tumor.obs["infiltration_type"].isin(["T", "B", "M"]), :
        ],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "infiltration_type",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Sum",
    },
    "patient_immune_infiltration": {
        "dataset": adata_primary_tumor,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "immune_infiltration",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Sum",
        "tools": ["dorothea", "progeny", "cytosig", "cpdb"],
    },
    "patient_immune_infiltration_treatment_coding": {
        "dataset": adata_primary_tumor,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "immune_infiltration",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Treatment('-')",
        "tools": ["dorothea", "progeny", "cytosig", "cpdb"],
    },
    "luad_lscc": {
        "dataset": adata_primary_tumor[
            adata_primary_tumor.obs["condition"].isin(["LUAD", "LSCC"]), :
        ],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "condition",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Treatment('LUAD')",
        "tools": ["dorothea", "progeny", "cytosig", "cpdb"],
    },
    "early_advanced": {
        "dataset": adata_primary_tumor,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "tumor_stage",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Treatment('early')",
    },
    "cell_types": {
        "dataset": adata_primary_tumor,
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "cell_type_major",
        "lm_covariate_str": "+ patient",
        "contrasts": "Sum",
        "cell_type_column": "ct_all",
        # need to manually split this up by cell-type, because here, cell-type is not a variable, but a condition.
        "tools": ["progeny", "cytosig"],
    },
}

# %%
comparisons = {
    k: v
    for k, v in comparisons.items()
    if k == "patient_immune_infiltration_treatment_coding"
}

# %% [markdown]
# ---

# %%
TOOLS = {
    "progeny": run_progeny,
    "dorothea": run_dorothea,
    "cytosig": run_cytosig,
}


def _prepare_subset(adata, tools):
    res = {}
    for tool in tools:
        res[tool] = TOOLS[tool](adata)
    return res


def prepare_dataset(
    id_,
    *,
    dataset,
    cell_type_column,
    tools=("progeny", "cytosig"),
    column_to_test,
    **kwargs,
):
    """Split anndata by cell-type and run the different signature enrichment
    methods"""
    tools = [t for t in tools if t in TOOLS]
    print(f"Preparing dataset for {id_}:")
    dataset = dataset[
        (dataset.obs[cell_type_column] != "other")
        & ~pd.isnull(dataset.obs[column_to_test])
        & (dataset.obs[column_to_test] != "nan"),
        :,
    ].copy()
    dataset.obs[column_to_test] = pd.Categorical(
        dataset.obs[column_to_test].astype(str)
    )
    print(f"\tSplitting anndata by {cell_type_column}:")
    adata_by_cell_type = sh.util.split_anndata(dataset, cell_type_column)
    res_by_cell_type = process_map(
        _prepare_subset, adata_by_cell_type.values(), itertools.repeat(tools)
    )
    all_adatas = {}
    for tool in tools:
        all_adatas[tool] = {}
        for ct, tmp_res in zip(adata_by_cell_type, res_by_cell_type):
            all_adatas[tool][ct] = tmp_res[tool]

    return all_adatas


def compare_signatures(
    id_,
    all_adatas,
    *,
    pseudobulk_group_by,
    column_to_test,
    lm_covariate_str,
    contrasts,
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
            contrasts=contrasts,
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
    contrasts,
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
    tmp_ad.obs.set_index("sample", inplace=True)
    tmp_ad.obs[column_to_test] = (
        dataset.obs.loc[:, ["sample", column_to_test]]
        .drop_duplicates()
        .set_index("sample")[column_to_test]
    )
    tmp_res = (
        sh.lm.test_lm(
            tmp_ad,
            f"~ C({column_to_test}, {contrasts}) {lm_covariate_str}",
            column_to_test,
            progress=True,
            contrasts=contrasts,
            n_jobs=48,
            chunksize=1200,
        )
        .set_index("variable")
        .join(ad_cpdb.var, how="left")
        .sort_values("pvalue")
        .dropna(how="any")
        .pipe(sh.util.fdr_correction)
        .pipe(sh.util.log2_fc)
    )
    tmp_res.to_csv(f"{artifact_dir}/differential_signature_{id_}_cpdb.tsv", sep="\t")
    return tmp_res


# %%
# sh.lm.test_lm(
#     datasets["patient_immune_infiltration_treatment_coding"]["cytosig"]["Tumor cells"],
#     "~C(immune_infiltration, Treatment('-')) + dataset",
#     "immune_infiltration",
#     contrasts="Treatment('-')",
# )

# %%
datasets = {}
for id_, config in comparisons.items():
    datasets[id_] = prepare_dataset(id_, **config)

# %%
results = {}
for id_, config in comparisons.items():
    results[id_] = compare_signatures(id_, datasets[id_], **config)

# %%
for id_, config in comparisons.items():
    if "cpdb" in config.get("tools", ["cpdb"]):
        print(id_)
        results[id_]["cpdb"] = compare_cpdb(id_, **config)

# %%
# id_ = "luad_lscc"
# results[id_] = compare_signatures(id_, datasets[id_], **comparisons[id_])

# %%
# results[id_]["progeny"]

# %% [markdown]
# ### test cpdb cell-type vs cell-type

# %%
tmp_ad = ad_cpdb[
    ad_cpdb.obs["sample"].isin(adata_primary_tumor.obs["sample"].unique().tolist()), :
]
tmp_ad2 = sc.AnnData(
    X=tmp_ad.var.join(
        pd.DataFrame(tmp_ad.X.T, columns=tmp_ad.obs_names, index=tmp_ad.var_names)
    )
    .pivot(index=["source", "target", "cluster_2"], columns=["cluster_1"])
    .T
)
tmp_ad2.obs.reset_index(inplace=True)
tmp_ad2.var.reset_index(inplace=True)
tmp_ad2.var.index = ["_".join([ct, s, t]) for _, s, t, ct in tmp_ad2.var.itertuples()]
tmp_ad2.obs.columns = ["sample", "cell_type_major"]

# %%
tmp_res = (
    sh.lm.test_lm(
        tmp_ad2,
        "~ C(cell_type_major, Sum) + sample",
        "cell_type_major",
        progress=True,
        contrasts="Sum",
        n_jobs=48,
        chunksize=600,
    )
    .set_index("variable")
    .join(tmp_ad2.var, how="left")
    .sort_values("pvalue")
    .dropna(how="any")
    .pipe(sh.util.fdr_correction)
    .pipe(sh.util.log2_fc)
)
tmp_res.to_csv(f"{artifact_dir}/differential_signature_cell_types_cpdb.tsv", sep="\t")
results["cell_types"]["cpdb"] = tmp_res

# %% [markdown]
# ## Differentially expressed dorothea TFs
#
#  * nothing significant for Neutrophils

# %%
results["luad_lscc"]["dorothea"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(sh.lm.plot_lm_result_altair, title="TFs (tumor cells LUAD/LSCC)")

# %%
results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].to_csv("/home/sturm/Downloads/neutrophils_tfs.tsv", sep="\t")

# %%
results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].pipe(sh.util.fdr_correction).pipe(
    sh.lm.plot_lm_result_altair, title="TFs (Neutrophils tumor/normal)"
)

# %%
pb_dorothea = sh.pseudobulk.pseudobulk(
    datasets["tumor_normal"]["dorothea"]["Neutrophils"],
    groupby=["dataset", "patient", "origin"],
    aggr_fun=np.mean,
)

# %%
tfoi = results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
]["variable"][:30]

# %%
sh.pairwise.plot_paired_fc(
    pb_dorothea, groupby="origin", paired_by="patient", metric="diff", var_names=tfoi
).properties(height=150)

# %%
sh.pairwise.plot_paired(
    pb_dorothea, groupby="origin", paired_by="patient", var_names=tfoi
)

# %% [markdown]
# ## Differentially expressed progeny pathways in tumor cells

# %%
results["patient_immune_infiltration_treatment_coding"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    sh.lm.plot_lm_result_altair, title="Differential pathways (tumor cells)"
)

# %%
results["luad_lscc"]["progeny"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(sh.lm.plot_lm_result_altair, title="Differential pathways (tumor cells)")

# %% [markdown]
# ## Differential cytokine signalling in selected cell-types

# %%
results["patient_immune_infiltration_treatment_coding"]["cytosig"]["group"].unique()

# %%
tmp_cytosig = (
    results["patient_immune_infiltration_treatment_coding"]["cytosig"]
    .loc[
        lambda x: x["cell_type"].isin(["Tumor cells", "Stromal", "Endothelial cell"]), :
    ]
    .pipe(sh.util.fdr_correction)
)

# %%
for ct in tmp_cytosig["cell_type"].unique():
    try:
        sh.lm.plot_lm_result_altair(
            tmp_cytosig.loc[lambda x: x["cell_type"] == ct], title=f"Cytosig for {ct}"
        ).display()
    except AttributeError:
        pass

# %%
results["luad_lscc"]["cytosig"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(sh.lm.plot_lm_result_altair, title="Cytosig (tumor cells)")

# %% [markdown]
# ---

# %%
marker_genes = {}
for tf in regulons.columns:
    marker_genes[tf] = regulons[tf][lambda x: x != 0]

# %%
pd.concat(marker_genes).reset_index(name="direction").rename(
    columns={"level_0": "TF"}
).to_csv("/home/sturm/Downloads/dorothea_signature.csv")

# %%
marker_genes["SOX2"].index.tolist()

# %%
regulons["SOX2"][lambda x: x != 0]

# %%

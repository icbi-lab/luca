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
import infercnvpy as cnv
import pyreadr

alt.data_transformers.disable_max_rows()

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=anndata.ImplicitModificationWarning)

# %%
ah = AnnotationHelper()

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
patient_strat = pd.read_csv(
    nxfvars.get(
        "stratification_csv",
        "../../data/30_downstream_analyses/stratify_patients/stratification/artifacts/patient_stratification.csv",
    ),
    index_col=0,
)
split_anndata_dir = Path(
    nxfvars.get(
        "split_anndata_dir", "../../data/30_downstream_analyses/infercnv/split_anndata/"
    )
)
scevan_dir = Path(
    nxfvars.get("scevan_dir", "../../data/30_downstream_analyses/infercnv/scevan/")
)
MIN_TUMOR_CELLS = 10
cpus = nxfvars.get("cpus", 16)

# %%
threadpool_limits(cpus)

# %%
adata = sc.read_h5ad(path_adata)

# %% [markdown]
# ## Read SCEVAN results

# %%
scevan_dirs = [
    x.name
    for x in scevan_dir.glob("*")
    if x.is_dir() and x.name.startswith("full_atlas_merged")
]


# %%
def read_scevan(scevan_id):
    tmp_ad = sc.read_h5ad(split_anndata_dir / f"{scevan_id}.h5ad")
    scevan_res_table = scevan_dir / scevan_id / "scevan_result.csv"

    if not scevan_res_table.exists():
        print(f"It seems sceavan was not successful for {scevan_id}. Skipping!")
        return None

    anno = pyreadr.read_r(
        scevan_dir / scevan_id / f"{scevan_id}_count_mtx_annot.RData"
    )["count_mtx_annot"]

    try:
        cnv.io.read_scevan(tmp_ad, scevan_dir / scevan_id, scevan_res_table)
    except KeyError as e:
        if not "CNA_mtx_relat" in str(e):
            raise
        else:
            print(f"Could not load CNA_mtx for {scevan_id}. Skipping!")
            return None

    return sc.AnnData(X=tmp_ad.obsm["X_scevan"], var=anno, obs=tmp_ad.obs)


# %%
cnv_adatas = process_map(read_scevan, scevan_dirs, max_workers=cpus)

# %% [markdown]
# ## Create combined CNV anndata object

# %%
cnv_adatas = [x for x in cnv_adatas if x is not None]

# %%
cnv_ad = anndata.concat(cnv_adatas, join="outer", fill_value=0)

# %%
# no idea why anndata.merge does not keep var columns
unique_var = pd.concat([x.var for x in cnv_adatas]).drop_duplicates()

# %%
assert set(cnv_ad.var_names) == set(unique_var.index)

# %%
assert cnv_ad.obs_names.is_unique
assert cnv_ad.var_names.is_unique

# %%
cnv_ad.var = unique_var.reindex(cnv_ad.var_names)

# %%
# sort by position
cnv_ad = cnv_ad[:, cnv_ad.var.sort_values(["seqnames", "start"]).index].copy()

# %%
cnv_ad.obs["immune_infiltration"] = cnv_ad.obs.join(
    patient_strat.loc[:, ["patient", "immune_infiltration"]].set_index("patient"),
    on="patient",
)["immune_infiltration"]

# %%
# nans are ok, as not every patient has a stratum assigned
cnv_ad.obs["immune_infiltration"].value_counts(dropna=False)

# %%
tumor_cells_per_patient = (
    cnv_ad.obs.loc[lambda x: x["cell_type_major"] == "Tumor cells"]
    .groupby("patient")
    .size()
)

# %%
selected_patients = tumor_cells_per_patient[
    tumor_cells_per_patient > MIN_TUMOR_CELLS
].index

# %%
# consensus method: only keep "buckets" that occur in at least 75% of patients
all_var = pd.concat([x.var.assign(patient=x.obs["patient"][0]) for x in cnv_adatas])

# %%
patients_per_segment = (
    all_var.loc[:, ["gene_id", "patient"]].drop_duplicates().groupby("gene_id").size()
)
keep_var = patients_per_segment[
    patients_per_segment > 0.75 * len(selected_patients)
].index

# %%
cnv_ad_tumor = cnv_ad[
    ~cnv_ad.obs["immune_infiltration"].isnull()
    & (cnv_ad.obs["cell_type_major"] == "Tumor cells")
    & (cnv_ad.obs["patient"].isin(selected_patients)),
    keep_var,
].copy()

# %%
cnv_pseudobulk = sh.pseudobulk.pseudobulk(
    cnv_ad_tumor,
    groupby=["dataset", "condition", "tumor_stage", "immune_infiltration", "patient"],
    aggr_fun=np.mean,
    min_obs=0,  # already filtered in previous step ("selected_patients")
)

# %%
cnv_pseudobulk

# %% [markdown]
# ## Run comparison

# %%
means_per_group = sh.pseudobulk.pseudobulk(
    cnv_pseudobulk, groupby="immune_infiltration", aggr_fun=np.mean, min_obs=0
)

# %%
means_per_group.obs.set_index("immune_infiltration", inplace=True)

# %%
# find segments that have means in opposing directions in at least one group
segments_opposing_directions = (np.sum(means_per_group.X > 0, axis=0) > 0) & (
    np.sum(means_per_group.X < 0, axis=0) > 0
)

# %%
segments_passing_thres =     np.max(np.abs(means_per_group.X), axis=0) > 0.01

# %%
with warnings.catch_warnings():
    warnings.filterwarnings(action="ignore")
    res_robust = sh.compare_groups.lm.test_lm(
        cnv_pseudobulk[
            cnv_pseudobulk.obs["condition"].isin(["LUAD", "LUSC"]),
            segments_passing_thres & segments_opposing_directions,
        ].copy(),
        groupby="immune_infiltration",
        formula="~ C(immune_infiltration, Sum) + dataset + condition + tumor_stage",
        contrasts="Sum",
        robust=True,
        n_jobs=cpus,
    )

# %%
segmeans = means_per_group.to_df().T
segmeans.columns = [f"segmean_{x}" for x in segmeans.columns]

# %%
res_robust = res_robust.merge(
    cnv_pseudobulk.var, how="inner", left_on="variable", right_index=True
).merge(segmeans, how="left", left_on="variable", right_index=True)

# %%
res_robust.sort_values("coef")

# %%
res_robust.to_csv("/home/sturm/Downloads/cnv_lm_res_unfiltered.csv")

# %%
res_co1 = (
    res_robust.loc[lambda x: np.abs(x["coef"]) > 0.01]
    .pipe(sh.util.fdr_correction)
    .loc[lambda x: x["fdr"] < 0.1]
    .sort_values("pvalue")
)
res_co1.to_csv(
    "/home/sturm/Downloads/cnv_lm_res_filtered_abs_diff_gt_0.01_fdr_lt_0.1.csv"
)
res_co1

# %%
res_co2 = (
    res_robust.loc[lambda x: np.abs(x["coef"]) > 0.05]
    .pipe(sh.util.fdr_correction)
    .loc[lambda x: x["fdr"] < 0.1]
    .sort_values("pvalue")
)
res_co2.to_csv(
    "/home/sturm/Downloads/cnv_lm_res_filtered_abs_diff_gt_0.05_fdr_lt_0.1.csv"
)
res_co2

# %% [markdown]
# ## cellphonedb genes

# %%
cpdb = pd.read_csv("../../tables/cellphonedb_2022-04-06.tsv", sep="\t")

# %%
res_co2.merge(cpdb, how="inner", left_on="gene_name", right_on="source_genesymbol")

# %%
res_co2.merge(cpdb, how="inner", left_on="gene_name", right_on="target_genesymbol")

# %%

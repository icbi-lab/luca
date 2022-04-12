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

alt.data_transformers.disable_max_rows()

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=anndata.ImplicitModificationWarning)

# %%
threadpool_limits(16)

# %%
ah = AnnotationHelper()

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)

# %%
patient_strat = pd.read_csv(
    nxfvars.get(
        "stratification_csv",
        "../../data/30_downstream_analyses/stratify_patients/artifacts/patient_stratification.csv",
    ),
    index_col=0,
)
MIN_TUMOR_CELLS = 50

# %%
adata = sc.read_h5ad(path_adata)

# %% [markdown]
# # CNV

# %%
patient_strat = patient_strat.assign(
    patient_lc=lambda x: x["patient"].str.lower()
).set_index("patient_lc")

# %%
obs_df

# %%
dfs = {}
obs_list = []
for group in ["ID", "M", "mixed", "T"]:
    dfs[group] = pd.read_csv(
        f"/home/sturm/Downloads/cnv/spread_{group}.tsv", sep="\t"
    ).set_index(["chromosome", "start", "end"])
    obs_list.extend(
        [
            [
                p,
                group,
                p.replace("Immune_desert_", "")
                .replace("T_", "")
                .replace("M_", "")
                .replace("mixed_", ""),
            ]
            for p in dfs[group].columns
        ]
    )

# %%
cnv_df = (
    dfs["ID"]
    .join(dfs["M"], how="outer")
    .join(dfs["mixed"], how="outer")
    .join(dfs["T"], how="outer")
    .fillna(0)
)

# %%
cnv_ad = sc.AnnData(cnv_df.T)

# %%
obs_df = pd.DataFrame(
    obs_list, columns=["patient", "_group_george", "patient_lc"]
).set_index("patient")

# %%
cnv_ad.obs["patient"] = obs_df["patient_lc"]
cnv_ad.obs.set_index("patient", inplace=True)

# %%
cnv_ad.obs["dataset"] = patient_strat["dataset"]
cnv_ad.obs["group"] = patient_strat["immune_infiltration"]
cnv_ad.obs["condition"] = patient_strat["tumor_type_annotated"]
cnv_ad.obs["tumor_stage"] = patient_strat["tumor_stage"]

# %%
cnv_ad.obs.loc[lambda x: x["dataset"].isnull()]

# %%
cnv_ad.obs

# %%
cnv_ad.var = (
    cnv_ad.var.reset_index()
    .assign(idx=lambda x: x["chromosome"].astype(str) + "_" + x["start"].astype(str))
    .set_index("idx")
)

# %%
cnv_ad.var

# %%
means_per_group = sh.pseudobulk.pseudobulk(
    cnv_ad, groupby="group", aggr_fun=np.mean, min_obs=0
)

# %%
means_per_group.obs.set_index("group", inplace=True)

# %%
segments_passing_thres = means_per_group.var.index[
    np.max(np.abs(means_per_group.X), axis=0) > 0.05
]

# %%
res_robust = sh.compare_groups.lm.test_lm(
    cnv_ad[cnv_ad.obs["condition"].isin(["LUAD", "LUSC"]), :].copy(),
    groupby="group",
    formula="~ C(group, Treatment('desert')) + dataset + condition + tumor_stage",
    contrasts="Treatment('desert')",
    robust=True,
)

# %%
res_robust

# %%
segmeans = means_per_group.to_df().T
segmeans.columns = [f"segmean_{x}" for x in segmeans.columns]

# %%
res_robust = res_robust.merge(
    cnv_ad.var.reset_index(), how="inner", left_on="variable", right_on="idx"
).merge(segmeans.reset_index(), how="left", left_on="variable", right_on="idx")

# %%
res_robust

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

# %%

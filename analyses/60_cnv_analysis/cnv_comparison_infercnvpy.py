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
scevan_res_dirs = list(
    Path("../../data/30_downstream_analyses/infercnv/scevan/").glob(
        "**/output/*CNAmtx.RData"
    )
)
infercnvpy_res_files = list(
    Path("../../data/30_downstream_analyses/infercnv/infercnvpy/").glob("**/*.h5ad")
)

# %%
adata.obs["patient_lower"] = adata.obs["patient"].str.lower()

# %%
adatas_cnv = sh.util.split_anndata(adata, "patient_lower")

# %%
adatas_infercnvpy = {}
for infercnvpy_path in tqdm(infercnvpy_res_files):
    patient = str(infercnvpy_path).split("/")[-2].replace("full_atlas_merged_", "")
    adatas_infercnvpy[patient] = sc.read_h5ad(infercnvpy_path)

# %% [markdown]
# ## Test based on infercnvpy

# %%
patient_strat = patient_strat.assign(
    patient_lc=lambda x: x["patient"].str.lower()
).set_index("patient_lc")


# %%
def _get_cnv(adata):
    if np.sum(adata.obs["cell_type_major"] == "Tumor cells") < 50:
        return np.full((adata.obsm["X_cnv"].shape[1],), np.nan)
    else:
        return np.mean(
            adata[adata.obs["cell_type_major"] == "Tumor cells", :].obsm["X_cnv"],
            axis=0,
        ).A1


# %%
cnv_df = pd.DataFrame.from_dict(
    {k: _get_cnv(adata) for k, adata in adatas_infercnvpy.items()}, orient="index"
)

# %%
cnv_df = cnv_df.dropna(how="any")

# %%
cnv_ad = sc.AnnData(cnv_df)

# %%
cnv_ad.obs["dataset"] = patient_strat["dataset"]
cnv_ad.obs["group"] = patient_strat["immune_infiltration"]
cnv_ad.obs["condition"] = patient_strat["tumor_type_annotated"]
cnv_ad.obs["tumor_stage"] = patient_strat["tumor_stage"]

# %%
cnv_ad.obs

# %%
res_robust = sh.compare_groups.lm.test_lm(
    cnv_ad[cnv_ad.obs["condition"].isin(["LUAD", "LUSC"]), :].copy(),
    groupby="group",
    formula="~ C(group, Treatment('desert')) + dataset + condition + tumor_stage",
    contrasts="Treatment('desert')",
    robust=False,
)

# %%
res_robust = sh.compare_groups.lm.test_lm(
    cnv_ad[cnv_ad.obs["condition"].isin(["LUAD", "LUSC"]), :].copy(),
    groupby="group",
    formula="~ C(group, Sum) + dataset + condition + tumor_stage",
    contrasts="Sum",
    robust=True,
)

# %%
res_robust.sort_values("pvalue")

# %%
pd.set_option("display.max_rows", 300)

# %%
res_co1 = (
    res_robust.loc[lambda x: np.abs(x["coef"]) > 0.05]
    .pipe(sh.util.fdr_correction)
    .loc[lambda x: x["fdr"] < 0.1]
    .sort_values("pvalue")
)
res_co1

# %%

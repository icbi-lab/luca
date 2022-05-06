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

# %%
for scevan_dir in tqdm(scevan_res_dirs):
    patient = str(scevan_dir).split("/")[-3].replace("full_atlas_merged_", "")
    try:
        cnv.io.read_scevan(
            adatas_cnv[patient],
            scevan_dir.parent,
            subclones=False,
            scevan_res_table=scevan_dir.parent.parent / "scevan_result.csv",
        )
    except (ValueError, KeyError) as e:
        warnings.warn(f"Patient {patient} failed: " + str(e))


# %%
ithcna_infercnvpy = {}
for patient, tmp_ad in tqdm(adatas_infercnvpy.items()):
    ithcna_infercnvpy[patient] = cnv.tl.ithcna(
        tmp_ad, groupby="cell_type_major", inplace=False
    )

# %%
ithgex = {}
ithcna = {}
n_cells = {}
for patient, tmp_ad in tqdm(adatas_cnv.items()):
    if "X_scevan" not in tmp_ad.obsm:
        continue

    tmp_ad = tmp_ad[tmp_ad.obs["origin"] == "tumor_primary", :]
    ithgex[patient] = cnv.tl.ithgex(tmp_ad, groupby="cell_type_major", inplace=False)
    ithcna[patient] = cnv.tl.ithcna(
        tmp_ad, groupby="cell_type_major", inplace=False, use_rep="X_scevan"
    )
    n_cells[patient] = tmp_ad.obs["cell_type_major"].value_counts().to_dict()

# %%
diversity_per_patient = (
    (
        pd.DataFrame.from_dict(ithgex)
        .T.reset_index()
        .melt(id_vars="index", var_name="cell_type_major", value_name="ITHGEX")
        .merge(
            pd.DataFrame.from_dict(ithcna)
            .T.reset_index()
            .melt(id_vars="index", var_name="cell_type_major", value_name="ITHCNA"),
            on=["index", "cell_type_major"],
            how="outer",
        )
        .merge(
            pd.DataFrame.from_dict(n_cells)
            .T.reset_index()
            .melt(
                id_vars="index", var_name="cell_type_major", value_name="n_tumor_cells"
            ),
            on=["index", "cell_type_major"],
            how="outer",
        )
        .merge(
            pd.DataFrame.from_dict(ithcna_infercnvpy)
            .T.reset_index()
            .melt(
                id_vars="index",
                var_name="cell_type_major",
                value_name="ithcna_infercnvpy",
            ),
            on=["index", "cell_type_major"],
            how="outer",
        )
    )
    .loc[lambda x: x["cell_type_major"] == "Tumor cells", :]
    .set_index("index")
    .drop(columns="cell_type_major")
)

# %%
diversity_per_patient

# %%
patient_strat["patient_lc"] = patient_strat["patient"].str.lower()
patient_strat.set_index("patient_lc", inplace=True)

# %%
patient_strat2 = patient_strat.join(
    diversity_per_patient.loc[lambda x: x["n_tumor_cells"] >= MIN_TUMOR_CELLS, :],
    how="inner",
)
patient_strat2["dataset"] = patient_strat2["dataset"].astype(str)

# %%
ax = sns.scatterplot(
    x="ITHCNA", y="ithcna_infercnvpy", data=patient_strat2, hue="study"
)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

# %%
mod = smf.ols(
    "ITHCNA ~ ithcna_infercnvpy + dataset",
    data=patient_strat2,
)
res = mod.fit()
res.summary()

# %%
for var in ["ITHCNA", "ithcna_infercnvpy", "ITHGEX"]:
    # for var in ["cnvsum"]:
    print(var)
    fig, (ax3, ax4, ax5) = plt.subplots(
        1, 3, figsize=(15, 5), gridspec_kw={"wspace": 0.5}
    )

    sns.swarmplot(
        x="tumor_type_annotated", y=var, data=patient_strat2, ax=ax3, hue="study"
    )
    sns.boxplot(x="tumor_type_annotated", y=var, data=patient_strat2, ax=ax3, width=0.2)
    ax3.get_legend().remove()

    sns.swarmplot(
        x="immune_infiltration", y=var, data=patient_strat2, ax=ax4, hue="study"
    )
    sns.boxplot(x="immune_infiltration", y=var, data=patient_strat2, ax=ax4, width=0.8)
    ax4.get_legend().remove()

    sns.swarmplot(x="tumor_stage", y=var, data=patient_strat2, ax=ax5, hue="study")
    sns.boxplot(x="tumor_stage", y=var, data=patient_strat2, ax=ax5, width=0.8)
    ax5.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    # ax5.get_legend().remove()
    plt.show()
    
    fig.savefig(f"{artifact_dir}/ith_analysis_{var}.pdf", bbox_inches="tight")

# %%
mod = smf.ols(
    "ITHCNA ~ C(immune_infiltration) + tumor_type_annotated + tumor_stage + dataset + n_tumor_cells",
    data=patient_strat2,
)
mod.fit().summary()

# %%
mod = smf.ols(
    "ITHGEX ~ C(immune_infiltration) + tumor_type_annotated + tumor_stage + dataset + n_tumor_cells",
    data=patient_strat2,
)
mod.fit().summary()

# %% [markdown]
# # Correlation with cell-types

# %%
# only on primary tumor samples;
# exclude datasets with only a single cell-type
frac_by_condition = (
    adata.obs.loc[
        lambda x: (x["origin"] == "tumor_primary")
        & ~x["dataset"].isin(["Guo_Zhang_2018"])
    ]
    .groupby(["patient"], observed=True)
    .apply(lambda x: x.value_counts("cell_type_major", normalize=True))
    .reset_index(name="frac_cells")
)

# %%
# only on primary tumor samples;
# exclude datasets with only a single cell-type
cells_per_dataset = (
    adata.obs.loc[
        lambda x: (x["origin"] == "tumor_primary")
        & ~x["dataset"].isin(["Guo_Zhang_2018"])
    ]
    .groupby(["dataset"], observed=True)
    .apply(lambda x: x.value_counts("cell_type_major", normalize=False))
    .reset_index(name="n_cells")
)

# %%
cna_ct_df = patient_strat2.merge(frac_by_condition, on="patient")

# %%
cells_per_dataset

# %%
results = []
for x in ["ITHCNA", "ITHGEX", "ithcna_infercnvpy"]:
    for ct in sorted(cna_ct_df["cell_type_major"].unique()):
        tmp_data = cna_ct_df.loc[lambda x: x["cell_type_major"] == ct]
        # only keep datasets that have at least 100 cells of a type
        tmp_data = tmp_data.loc[
            lambda x: x["dataset"].isin(
                cells_per_dataset.loc[
                    lambda x: (x["n_cells"] > 50) & (x["cell_type_major"] == ct),
                    "dataset",
                ]
            ),
            [x, "frac_cells", "dataset", "n_tumor_cells"],
        ].dropna(how="any")
        fig, ax = plt.subplots(1, 1)
        ax = sns.scatterplot(
            data=tmp_data,
            x=x,
            y="frac_cells",
            hue="dataset",
        )
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        ax.set_ylabel(ct)
        plt.show()
        res = smf.rlm(
            f"{x} ~ frac_cells + dataset + n_tumor_cells", data=tmp_data
        ).fit()
        results.append(
            [
                x,
                ct,
                res.params["Intercept"],
                res.params["frac_cells"],
                res.pvalues["frac_cells"],
                *scipy.stats.pearsonr(tmp_data[x], tmp_data["frac_cells"]),
            ]
        )

# %%
res_df = (
    pd.DataFrame.from_records(
        results,
        columns=[
            "metric",
            "cell_type",
            "intercept",
            "coef",
            "pvalue",
            "pearson",
            "pvalue_pearson",
        ],
    )
    .pipe(sh.util.fdr_correction)
    .sort_values("fdr")
)

# %%
res_df[:20]

# %%
res_df.pipe(
    sh.compare_groups.pl.plot_lm_result_altair,
    x="cell_type",
    y="metric",
    color="pearson",
    p_cutoff=1,
    reverse=True,
)

# %%

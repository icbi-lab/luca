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
    "../../data/30_downstream_analyses/04_neutrophil_subclustering//artifacts/full_atlas_neutrophil_clusters.h5ad",
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
# For the patient stratification, treat the two batches of the UKIM-V dataset as one
adata.obs["dataset"] = adata.obs["dataset"].str.replace("UKIM-V-2", "UKIM-V")

# %%
patient_strat.set_index("patient", inplace=True)
patient_strat["dataset"] = (
    adata.obs.loc[:, ["patient", "dataset"]]
    .drop_duplicates()
    .set_index("patient")["dataset"]
)
patient_strat.reset_index(inplace=True)

# %% [markdown]
# # CNV

# %%
# TODO: compare to SCEVAN
# TODO: regress out dataset-specific effects, or at least include dataset in linear model.

# %%
# scevan_res_files = list(
#     Path("../../data/30_downstream_analyses/infercnv/scevan/").glob(
#         "**/scevan_result.csv"
#     )
# )
adatas_infercnvpy = list(
    Path("../../data/30_downstream_analyses/infercnv/infercnvpy/").glob("**/*.h5ad")
)

# %%
adatas_cnv = {}
for f in tqdm(adatas_infercnvpy):
    patient = str(f).split("/")[-2].replace("full_atlas_merged_", "")
    adatas_cnv[patient] = sc.read_h5ad(f)

# %%
ithgex = {}
ithcna = {}
for patient, tmp_ad in tqdm(adatas_cnv.items()):
    ithgex[patient] = sh.diversity.ithgex(
        tmp_ad, groupby="cell_type_major", inplace=False
    )
    ithcna[patient] = sh.diversity.ithcna(
        tmp_ad, groupby="cell_type_major", inplace=False
    )

# %%
cnv_scores = {}
for patient, tmp_ad in tqdm(adatas_cnv.items()):
    cnv_scores[patient] = sh.diversity.cnv_score(
        tmp_ad, obs_key="cell_type_major", inplace=False
    )

# %%
n_tumor_cells = {}
for patient, tmp_ad in tqdm(adatas_cnv.items()):
    n_tumor_cells[patient] = tmp_ad.obs["cell_type_major"].value_counts().to_dict()

# %%
diversity_per_patient = (
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
        pd.DataFrame.from_dict(cnv_scores)
        .T.reset_index()
        .melt(id_vars="index", var_name="cell_type_major", value_name="cnv_score"),
        on=["index", "cell_type_major"],
        how="outer",
    )
    .merge(
        pd.DataFrame.from_dict(n_tumor_cells)
        .T.reset_index()
        .melt(id_vars="index", var_name="cell_type_major", value_name="n_tumor_cells"),
        on=["index", "cell_type_major"],
        how="outer",
    )
)

# %%
ith_df = diversity_per_patient.loc[
    lambda x: x["cell_type_major"] == "Tumor cells", :
].set_index("index")

# %%
patient_strat["patient_lc"] = patient_strat["patient"].str.lower()
patient_strat.set_index("patient_lc", inplace=True)

# %%
ith_df.loc[lambda x: x["n_tumor_cells"] >= 100, :]

# %%
patient_strat2 = patient_strat.join(
    ith_df.loc[lambda x: x["n_tumor_cells"] >= 50, :], how="inner"
)
patient_strat2["dataset"] = patient_strat2["dataset"].astype(str)

# %%
# TODO use linear model to check for confounding effects: dataset and number of cells!
for var in ["ITHCNA", "ITHGEX", "cnv_score"]:
    # for var in ["cnvsum"]:
    print(var)
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, figsize=(22, 5))
    sns.swarmplot(x="TMIG", y=var, data=patient_strat2, ax=ax1, hue="dataset")
    sns.boxplot(x="TMIG", y=var, data=patient_strat2, ax=ax1, width=0.2)
    ax1.get_legend().remove()

    sns.swarmplot(
        x="infiltration_state", y=var, data=patient_strat2, ax=ax2, hue="dataset"
    )
    sns.boxplot(x="infiltration_state", y=var, data=patient_strat2, ax=ax2, width=0.2)
    ax2.get_legend().remove()

    sns.swarmplot(
        x="tumor_type_annotated", y=var, data=patient_strat2, ax=ax3, hue="dataset"
    )
    sns.boxplot(x="tumor_type_annotated", y=var, data=patient_strat2, ax=ax3, width=0.2)
    ax3.get_legend().remove()

    sns.swarmplot(
        x="immune_infiltration", y=var, data=patient_strat2, ax=ax4, hue="dataset"
    )
    sns.boxplot(x="immune_infiltration", y=var, data=patient_strat2, ax=ax4, width=0.8)
    ax4.get_legend().remove()

    sns.swarmplot(x="tumor_stage", y=var, data=patient_strat2, ax=ax5, hue="dataset")
    sns.boxplot(x="tumor_stage", y=var, data=patient_strat2, ax=ax5, width=0.8)
    ax5.get_legend().remove()
    plt.show()

# %%
mod = smf.ols(
    "ITHCNA ~ C(immune_infiltration) + dataset + n_tumor_cells", data=patient_strat2
)
mod.fit().summary()

# %%
mod = smf.ols(
    "ITHGEX ~ C(immune_infiltration) + dataset + n_tumor_cells", data=patient_strat2
)
mod.fit().summary()

# %%
mod = smf.ols(
    "cnv_score ~ C(immune_infiltration) + dataset + n_tumor_cells", data=patient_strat2
)
mod.fit().summary()

# %%
mod = smf.ols(
    "cnv_score ~ C(tumor_type_annotated) + dataset + n_tumor_cells",
    data=patient_strat2.loc[lambda x: x["tumor_type_annotated"].isin(["LUAD", "LUSC"])],
)
mod.fit().summary()


# %% [markdown]
# # Correlation with cell-types

# %%
def get_cell_type_group(ct):
    for group, cts in cell_types.items():
        if ct in cts:
            return group
    return "other"


# %%
adata_primary_tumor = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"]),
    :,
]

# %%
ad_cts = sc.AnnData(
    X=(
        adata_primary_tumor.obs.groupby(
            ["dataset", "patient", "cell_type_major"], observed=True
        )
        .size()
        .reset_index(name="n")
        .pivot_table(
            values="n",
            columns="cell_type_major",
            index=["dataset", "patient"],
            fill_value=0,
        )
    )
)
ad_cts.obs = ad_cts.obs.reset_index().set_index("patient")

# %%
sc.pp.normalize_total(ad_cts, target_sum=1)

# %%
sc.pl.matrixplot(
    ad_cts,
    var_names=ad_cts.var_names,
    groupby="patient",
    dendrogram=False,
    swap_axes=True,
    cmap="viridis",
    # vmin=-0.25,
    # vmax=0.25
    # # vmin=0,
    # vmax=1,
    standard_scale="var",
)

# %%
ad_cts.obs["patient"] = ad_cts.obs.index.values

# %%
ad_cts.obs.index = ad_cts.obs.index.str.lower()

# %%
ct_df = pd.DataFrame(ad_cts.X, index=ad_cts.obs_names, columns=ad_cts.var_names).join(
    ad_cts.obs
)

# %%
ith_df_subset = ith_df.loc[lambda x: x["n_tumor_cells"] >= 100, :]

# %%
cna_ct_df = ith_df_subset.join(ct_df, how="inner")

# %%
cna_ct_df

# %%
for x in ["ITHCNA", "ITHGEX"]:
    for y in ["T cell CD8", "Neutrophils", "Macrophage"]:
        fig, ax = plt.subplots(1, 1)
        ax = sns.scatterplot(data=cna_ct_df, x=x, y=y, hue="dataset")
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.show()
        mod = smf.ols(f"{x} ~ Q('{y}') + dataset", data=cna_ct_df)
        display(mod.fit().summary())

# %%

# %%

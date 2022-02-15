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
import numpy as np
import scanpy_helpers as sh

sc.settings.set_figure_params(figsize=(5, 5))

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
adata_primary_tumor = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])
    & (adata.obs["cell_type_major"] != "other"),
    :,
]

# %%
sc.pl.umap(adata_primary_tumor, color="cell_type_major")

# %% [markdown]
# ### Subsample
#
# For cibersortx, we subsample the dataset, to balance the effect of platforms, datasets and patients, 
# and at the same time, reduce the dataset size, s.t. it can be managed by cibersortx. 

# %%
adata_sub = sh.deconvolution.balanced_subsample(
    adata_primary_tumor,
    cell_type_key="cell_type_major",
    patient_key="patient",
    n_each=5,
)

# %%
adata_sub.obs["cell_type_major"].value_counts()

# %%
adata_sub.obs.loc[:, "patient"].to_csv(
    "/home/sturm/Downloads/cibersortx/batch_ids.tsv", sep="\t"
)

# %%
cibersort_df = pd.DataFrame(
    adata_sub.X.todense().T,
    index=adata_sub.var_names,
    columns=adata_sub.obs["cell_type_major"],
).astype(int)

# %%
cibersort_df.to_csv(
    "/home/sturm/Downloads/cibersort_single_cell_matrix_5.txt", sep="\t"
)

# %% [markdown]
# ## marker gene analysis

# %%
adata_primary_tumor.X.data[:20]

# %%
pb_primary_tumor = sh.pseudobulk.pseudobulk(
    adata_primary_tumor, groupby=["patient", "cell_type_major"]
)

# %%
pb_primary_tumor.obs

# %%
sc.pp.normalize_total(pb_primary_tumor, target_sum=1e6)
sc.pp.log1p(pb_primary_tumor)

# %%
pb_primary_tumor.var["fold_change"] = sh.signatures.fold_change(
    pb_primary_tumor, pb_primary_tumor.obs["cell_type_major"] == "Neutrophils"
)

# %%
pb_primary_tumor.var["specific_fold_change"] = sh.signatures.specific_fold_change(
    pb_primary_tumor, obs_col="cell_type_major", positive_class="Neutrophils"
)

# %%
(pb_primary_tumor.obs["cell_type_major"] == "Neutrophils").value_counts()

# %%
pb_primary_tumor.var["roc_auc"] = sh.signatures.roc_auc(
    pb_primary_tumor, pb_primary_tumor.obs["cell_type_major"] == "Neutrophils"
)

# %%
pb_primary_tumor.var

# %%
signature_genes = (
    pb_primary_tumor.var.query("roc_auc >=0.97")
    .query("fold_change > 3")
    .query("specific_fold_change > 2")
)
signature_genes

# %% [markdown]
# ## Marker gene analysis 2 (including normal; TAN signature)

# %%
adata_tumor = adata[
    adata.obs["condition"].isin(["LUAD", "LSCC"])
    & adata.obs["origin"].isin(["normal_adjacent", "normal", "tumor_primary"]),
    :,
].copy()

# %%
adata_tumor.obs["cell_type_tan"] = adata_tumor.obs["cell_type_major"].astype(str)

# %%
adata_tumor.obs["origin"].unique()

# %%
adata_tumor.obs.loc[
    (adata_tumor.obs["origin"] == "tumor_primary")
    & (adata_tumor.obs["cell_type_tan"] == "Neutrophils"),
    "cell_type_tan",
] = "TAN"

# %%
pb_adata_tumor = sh.pseudobulk.pseudobulk(adata_tumor, groupby=["patient", "cell_type_tan"])

# %%
pb_adata_tumor

# %%
sc.pp.normalize_total(pb_adata_tumor, target_sum=1e6)
sc.pp.log1p(pb_adata_tumor)

# %%
pb_adata_tumor.var["fold_change"] = sh.signatures.fold_change(
    pb_adata_tumor, pb_adata_tumor.obs["cell_type_tan"] == "TAN"
)

# %%
pb_adata_tumor.var["specific_fold_change"] = sh.signatures.specific_fold_change(
    pb_adata_tumor, obs_col="cell_type_tan", positive_class="TAN"
)

# %%
pb_adata_tumor.var["roc_auc"] = sh.signatures.roc_auc(
    pb_adata_tumor, pb_adata_tumor.obs["cell_type_tan"] == "TAN"
)

# %%
signature_genes = (
    pb_adata_tumor.var.query("roc_auc >=0.97")
    .query("fold_change > 2")
    .query("specific_fold_change > 1")
)
signature_genes

# %%
signature_genes.index.tolist()

# %%

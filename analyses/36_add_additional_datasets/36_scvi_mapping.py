# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: SSH apollo-15 pircher-sc-integrate2
#     language: ''
#     name: rik_ssh_apollo_15_pircherscintegrate2
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
import scvi

sc.settings.set_figure_params(figsize=(4, 4))
import pandas as pd
from nxfvars import nxfvars
from scanpy_helpers.integration import (
    sanitize_adata,
    validate_adata,
    aggregate_duplicate_gene_symbols,
)
import scipy.sparse
import anndata
import sklearn
from threadpoolctl import threadpool_limits
import matplotlib.pyplot as plt
import scanpy_helpers as sh
from scanpy_helpers.annotation import AnnotationHelper
import numpy as np
from tqdm import tqdm

# %% [markdown]
# ## Get Input data

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
dataset_table = pd.read_csv(
    nxfvars.get("samplesheet", "../../tables/samplesheet_scrnaseq_preprocessing2.csv")
)
dataset_path = nxfvars.get(
    "dataset_path", "../../data/30_downstream_analyses/01_qc_and_filtering/"
)

# %%
threadpool_limits(nxfvars.get("cpus", 16))
reference_scanvi = sc.read_h5ad(
    nxfvars.get(
        "reference_scanvi_h5ad",
        "../../data/20_build_atlas/annotate_datasets/35_final_atlas/full_atlas_hvg_integrated_scvi_integrated_scanvi.h5ad",
    )
)
# this is needed as the 'reference' above only contains the HVG for SCVI
reference_atlas = sc.read_h5ad(
    nxfvars.get(
        "reference_atlas",
        "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
    )
)

# %%
vae_ref_scanvi = scvi.model.SCANVI.load(
    nxfvars.get(
        "reference_scanvi_model",
        "../../data/20_build_atlas/annotate_datasets/35_final_atlas/full_atlas_hvg_integrated_scvi_scanvi_model/",
    ),
    adata=reference_scanvi,
)

# %%
gene_symbol_df = pd.read_csv(
    nxfvars.get("gene_symbol_table", "../../tables/gene_symbol_dict.csv"),
    index_col=False,
)
gene_symbol_dict = {
    alias: symbol
    for alias, symbol in zip(gene_symbol_df["alias"], gene_symbol_df["gene_symbol"])
}

# %%
datasets = {
    dataset_id: sc.read_h5ad(
        f"{dataset_id}.qc.h5ad"
        if dataset_path == "."
        else f"{dataset_path}/{dataset_id}/{dataset_id}.qc.h5ad"
    )
    for dataset_id in tqdm(dataset_table["id"])
}

# %%
sc.pl.umap(reference_scanvi, color="cell_type")

# %% [markdown]
# ## Standardize query dataset
# metadata, gene symbols, etc.

# %%
datasets["UKIM-V-2"].obs["dataset"] = "UKIM-V-2"
datasets["UKIM-V-2"].obs["patient"] = [
    f"UKIM-V-2_{p}" for p in datasets["UKIM-V-2"].obs["patient"]
]
datasets["UKIM-V-2"].obs["tissue"] = "lung"
datasets["UKIM-V-2"].obs["tumor_stage"] = np.nan
datasets["UKIM-V-2"].obs.loc[
    datasets["UKIM-V-2"].obs["uicc_stage"].isin(["I", "II", "IA"]), "tumor_stage"
] = "early"
datasets["UKIM-V-2"].obs.loc[
    datasets["UKIM-V-2"]
    .obs["uicc_stage"]
    .isin(["III", "IV", "IIIA", "IIIB", "III or IV"]),
    "tumor_stage",
] = "late"
assert np.all(
    pd.isnull(datasets["UKIM-V-2"].obs["uicc_stage"])
    == pd.isnull(datasets["UKIM-V-2"].obs["tumor_stage"])
)

# %%
for dataset_id, adata in datasets.items():
    adata.obs["dataset"] = dataset_id
    print(f"Validating {dataset_id}")
    sanitize_adata(adata)
    validate_adata(adata)

# %% [markdown]
# ### Gene identifier remapping

# %%
for id_, query in datasets.items():
    # remap gene symbols
    query.var_names = [gene_symbol_dict.get(x, x) for x in query.var_names]
    query = aggregate_duplicate_gene_symbols(query)
    query.obs = query.obs.loc[
        :, [x for x in query.obs.columns if x in reference_atlas.obs.columns]
    ]
    datasets[id_] = query

# %% [markdown]
# ### aggregate duplicate gene symbols

# %%
for dataset_id, dataset in datasets.items():
    print(dataset_id)
    datasets[dataset_id] = aggregate_duplicate_gene_symbols(dataset)

# %% [markdown]
# ### Log-transform and normalize

# %%
for id_, query in datasets.items():
    query_raw = query.copy()
    sc.pp.normalize_total(query_raw)
    sc.pp.log1p(query_raw)
    query.raw = query_raw

# %%
# set variables for scvi
for query in datasets.values():
    query.obs["batch"] = query.obs["sample"]
    query.obs["cell_type"] = "unknown"

# %%
# keep all gene symbols, will reindex later to same genes as reference dataset
query_merged = sh.integration.merge_datasets(datasets.values(), symbol_in_n_datasets=1)

# %%
query_merged = query_merged[:, reference_atlas.var_names].copy()

# %%
query_merged_hvg = query_merged[:, reference_scanvi.var_names].copy()

# %% [markdown]
# ## Integrate using scANVI
# (using the scArches approach, but as implemented in the scvi package)

# %%
# the model version warning doesn't matter see
# https://github.com/scverse/scvi-tools/pull/1431
vae_q = scvi.model.SCANVI.load_query_data(query_merged_hvg, vae_ref_scanvi)

# %%
vae_q.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

# %%
new_samples = query_merged_hvg.obs["sample"].unique()

# %%
query_merged_hvg.obs["doublet_status"] = None

@sh.util.suppress_stdout
def predict_solo(batch):
    scvi.settings.verbosity = scvi.logging.CRITICAL
    tmp_solo = scvi.external.SOLO.from_scvi_model(
        vae_q, restrict_to_batch=s
    )
    tmp_solo.train(train_size=0.9)
    doublet_prediction = tmp_solo.predict(False)
    # for whatever reason, another "-0" is appended to the index. 
    doublet_prediction.index = doublet_prediction.index.str.replace("-0$", "", regex=True)
    scvi.settings.verbosity = scvi.logging.INFO
    return doublet_prediction

for s in tqdm(new_samples):
    doublet_prediction = predict_solo(s)
    query_merged_hvg.obs.loc[lambda x: x["sample"] == s, "doublet_status"] = doublet_prediction

# %%
query_merged_hvg.obs["doublet_status"].value_counts(dropna=False)

# %%
query_merged_hvg.obsm["X_scANVI"] = vae_q.get_latent_representation()
query_merged_hvg.obs["_predictions"] = vae_q.predict()


# %% [markdown]
# ## Visualize query dataset

# %%
sc.pp.neighbors(query_merged_hvg, use_rep="X_scANVI")
sc.tl.leiden(query_merged_hvg, key_added="_leiden")
sc.tl.umap(query_merged_hvg)

# %%
sc.pl.umap(query_merged_hvg, color=["_predictions", "_leiden", "doublet_status"], wspace=1)

# %%
ah = AnnotationHelper()

# %% tags=[]
ah.plot_dotplot(query_merged_hvg, groupby="_predictions")

# %% [markdown]
# ## Merge query and reference

# %%
reference_atlas.obsm["X_scANVI"] = reference_scanvi.obsm["X_scANVI"]

# %%
adata_full = anndata.concat(
    [reference_atlas, query_merged_hvg],
    merge="first",
    uns_merge="first",
    join="outer",
)
assert adata_full.var_names.is_unique
assert adata_full.obs_names.is_unique

# %%
# raw does not get concatenated... let's recompute it
adata_full_raw = adata_full.copy()
sc.pp.normalize_total(adata_full_raw)
sc.pp.log1p(adata_full_raw)
adata_full.raw = adata_full_raw

# %%
# Excluding Maier_Merad_2020 which is a subset of Leader_Merad_2021
adata_full = adata_full[adata_full.obs["dataset"] != "Maier_Merad_2020", :]
adata_full = adata_full[adata_full.obs["doublet_status"] != "doublet", :].copy()

# %%
sc.pp.neighbors(adata_full, use_rep="X_scANVI")


# %%
# initalize umap with the original coordinates. Missing values (=query) are initialized with random values.
init_pos_df = pd.DataFrame(reference_scanvi.obsm["X_umap"], index=reference_scanvi.obs_names).reindex(
    adata_full.obs_names
)
for col in init_pos_df.columns:
    na_mask = init_pos_df[col].isnull()
    init_pos_df.loc[na_mask, col] = np.random.uniform(
        np.min(init_pos_df[col]), np.max(init_pos_df[col]), size=np.sum(na_mask)
    )

# %%
sc.tl.umap(adata_full, init_pos=init_pos_df.values)

# %% [markdown]
# ## Automated cell-type annotation
#
# The annotation by scANVI doesn't look very good. I, therefore, implemented an
# approach based on nearest neighbor majority voting.

# %%
annot_cols = ["cell_type", "cell_type_major", "cell_type_tumor", "cell_type_coarse"]

# %%
for col in annot_cols:
    # I want the prediction in a separate column...
    sh.annotation.classify_cell_types_nearest_neighbors(
        adata_full,
        col,
        mask_reference=~adata_full.obs["dataset"].isin(datasets.keys()),
        mask_query=adata_full.obs["dataset"].isin(datasets.keys()),
        key_added=f"_{col}_predicted",
    )
    # ... but also merged back into the original column
    adata_full.obs.loc[
        ~adata_full.obs[f"_{col}_predicted"].isnull(), col
    ] = adata_full.obs[f"_{col}_predicted"]

# %%
list(datasets.keys())

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(adata_full, color="dataset", groups=list(datasets.keys()), size=4)

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    for col in annot_cols:
        sc.pl.umap(adata_full, color=[col, f"_{col}_predicted"], size=2)

# %% [markdown]
# ## Write output
#

# %%
adata_full.write_h5ad(f"{artifact_dir}/full_atlas_merged.h5ad")

# %%

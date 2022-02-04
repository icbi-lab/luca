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

# %% [markdown]
# ## Get Input data

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
threadpool_limits(nxfvars.get("cpus", 16))
query = sc.read_h5ad(
    nxfvars.get(
        "query",
        "../../data/30_downstream_analyses/01_qc_and_filtering/UKIM-V-2/UKIM-V-2.qc.h5ad",
    )
)
reference = sc.read_h5ad(
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
vae_ref = scvi.model.SCANVI.load(
    nxfvars.get(
        "reference_scanvi_model",
        "../../data/20_build_atlas/annotate_datasets/35_final_atlas/full_atlas_hvg_integrated_scvi_scanvi_model/",
    ),
    adata=reference,
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
sc.pl.umap(reference, color="cell_type")

# %% [markdown]
# ## Standardize query dataset
# metadata, gene symbols, etc.

# %%
query.obs["dataset"] = "UKIM-V-2"
query.obs["patient"] = [f"UKIM-V-2_{p}" for p in query.obs["patient"]]
query.obs["tissue"] = "lung"
query.obs["tumor_stage"] = np.nan
query.obs.loc[query.obs["uicc_stage"].isin(["I", "II", "IA"]), "tumor_stage"] = "early"
query.obs.loc[
    query.obs["uicc_stage"].isin(["III", "IV", "IIIA", "IIIB", "III or IV"]),
    "tumor_stage",
] = "late"
assert np.all(pd.isnull(query.obs["uicc_stage"]) == pd.isnull(query.obs["tumor_stage"]))

# %%
sanitize_adata(query)

# %%
validate_adata(query)

# %%
# remap gene symbols
query.var_names = [gene_symbol_dict.get(x, x) for x in query.var_names]
query = aggregate_duplicate_gene_symbols(query)

# %%
query.obs = query.obs.loc[
    :, [x for x in query.obs.columns if x in reference_atlas.obs.columns]
]

# %%
query_raw = query.copy()
sc.pp.normalize_total(query_raw)
sc.pp.log1p(query_raw)
query.raw = query_raw

# %%
# There are a few gene symbols in the 6000 variable genes of the reference that don't appear in the query:
len(set(reference.var_names) - set(query.var_names))

# %%
# Therefore, we re-index the query anndata object to fill those genes with zeros
query_hvg = sh.util.reindex_adata(query, reference.var_names)

# %% [markdown]
# ## Integrate using scANVI
# (using the scArches approach, but as implemented in the scvi package)

# %%
# set variables for scvi
query_hvg.obs["batch"] = query_hvg.obs["sample"]
query_hvg.obs["cell_type"] = "unknown"

# %%
vae_q = scvi.model.SCANVI.load_query_data(query_hvg, vae_ref)

# %%
vae_q.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

# %%
query.obsm["X_scANVI"] = vae_q.get_latent_representation()
query.obs["_predictions"] = vae_q.predict()


# %% [markdown]
# ## Visualize query dataset

# %%
sc.pp.neighbors(query, use_rep="X_scANVI")
sc.tl.leiden(query, key_added="_leiden")
sc.tl.umap(query)

# %%
sc.pl.umap(query, color=["_predictions", "_leiden"], wspace=1)

# %%
ah = AnnotationHelper()

# %% tags=[]
ah.plot_dotplot(query, groupby="_predictions")

# %% [markdown]
# ## Merge query and reference

# %%
reference_atlas.obsm["X_scANVI"] = reference.obsm["X_scANVI"]

# %%
adata_full = anndata.concat(
    [reference_atlas, sh.util.reindex_adata(query, reference_atlas.var_names)],
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
sc.pp.neighbors(adata_full, use_rep="X_scANVI")


# %%
# initalize umap with the original coordinates. Missing values (=query) are initialized with random values.
init_pos_df = pd.DataFrame(reference.obsm["X_umap"], index=reference.obs_names).reindex(
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
        mask_reference=adata_full.obs["dataset"] != "UKIM-V-2",
        mask_query=adata_full.obs["dataset"] == "UKIM-V-2",
        key_added=f"_{col}_predicted",
    )
    # ... but also merged back into the original column
    adata_full.obs.loc[
        ~adata_full.obs[f"_{col}_predicted"].isnull(), col
    ] = adata_full.obs[f"_{col}_predicted"]

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(adata_full, color="dataset", groups=["UKIM-V-2"], size=4)

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

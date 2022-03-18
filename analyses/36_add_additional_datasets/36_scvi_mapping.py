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
dataset_id = nxfvars.get("dataset_id", "UKIM-V-2")
dataset_path = nxfvars.get(
    "dataset_path", "../../data/30_downstream_analyses/01_qc_and_filtering/"
)

# %%
threadpool_limits(nxfvars.get("cpus", 2))
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
query = sc.read_h5ad(
    f"{dataset_id}.qc.h5ad"
    if dataset_path == "."
    else f"{dataset_path}/{dataset_id}/{dataset_id}.qc.h5ad"
)

# %%
sc.pl.umap(reference_scanvi, color="cell_type")

# %% [markdown]
# ## Standardize query dataset
# metadata, gene symbols, etc.

# %%
sh.integration.validate_adata(query, validate_obs=False)

# %% [markdown]
# ### Gene identifier remapping

# %%
# remap gene symbols
query.var_names = [gene_symbol_dict.get(x, x) for x in query.var_names]
query = aggregate_duplicate_gene_symbols(query)
query.obs = query.obs.loc[
    :, [x for x in query.obs.columns if x in reference_atlas.obs.columns]
]

# %% [markdown]
# ### aggregate duplicate gene symbols

# %%
query = aggregate_duplicate_gene_symbols(query)

# %%
# set variables for scvi
query.obs["batch"] = query.obs["sample"]
query.obs["cell_type"] = "unknown"

# %%
query = sh.util.reindex_adata(query, reference_atlas.var_names)

# %%
query_hvg = query[:, reference_scanvi.var_names].copy()

# %% [markdown]
# ## Integrate using scANVI
# (using the scArches approach, but as implemented in the scvi package)

# %%
# the model version warning doesn't matter see
# https://github.com/scverse/scvi-tools/pull/1431
vae_q = scvi.model.SCANVI.load_query_data(query_hvg, vae_ref_scanvi)

# %%
vae_q.train(
    max_epochs=100,
    plan_kwargs=dict(weight_decay=0.0),
    check_val_every_n_epoch=10,
)

# %%
new_samples = query_hvg.obs["sample"].unique()

# %%
query_hvg.obs["doublet_status"] = None

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
    query_hvg.obs.loc[lambda x: x["sample"] == s, "doublet_status"] = doublet_prediction

# %%
query.obs["doublet_status"] = query_hvg.obs["doublet_status"]
query.obs["doublet_status"].value_counts(dropna=False)

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
sc.pl.umap(query, color=["_predictions", "_leiden", "doublet_status"], wspace=1)

# %%
ah = AnnotationHelper()

# %% tags=[]
ah.plot_dotplot(query, groupby="_predictions")

# %% [markdown]
# # Write output

# %%
# makes issues while saving to anndata 
query.var.drop("mito", axis="columns", errors="ignore", inplace=True)

# %%
query.write_h5ad(f"{artifact_dir}/{dataset_id}_integrated.h5ad")

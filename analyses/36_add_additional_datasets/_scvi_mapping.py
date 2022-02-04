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
from scanpy_helpers.integration import sanitize_adata, validate_adata, aggregate_duplicate_gene_symbols
import scipy.sparse
import anndata

# %%
test = sc.read_h5ad("../../data/20_build_atlas/annotate_datasets/31_cell_types_coarse/artifacts/adata_cell_type_coarse.h5ad")

# %%
sc.pl.umap(test, color="cell_type")

# %%
query = sc.read_h5ad(
    "../../data/30_downstream_analyses/01_qc_and_filtering/UKIM-V-2/UKIM-V-2.qc.h5ad"
)

# %%
# reference = sc.read_h5ad("../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad")
reference = sc.read_h5ad(
    "../../data/20_build_atlas/integrate_datasets/24_scanvi_umap/all.umap_leiden.h5ad"
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

# %%
vae_ref = scvi.model.SCANVI.load(
    "../../data/20_build_atlas/integrate_datasets/23_scanvi/all_hvg_integrated_scvi_scanvi_model/",
    adata=reference,
)

# %%
query.var

# %%
query.obs["dataset"] = "UKIM-V-2"
query.obs["tissue"] = "lung"

# %%
sanitize_adata(query)

# %%
validate_adata(query)

# %%
query.var_names = [gene_symbol_dict.get(x, x) for x in query.var_names]

# %%
query = aggregate_duplicate_gene_symbols(query)

# %%
set(reference.var_names) - set(query.var_names)


# %%
def reindex_adata(adata, new_var_names):
    tmp_ad =  sc.AnnData(
            var = adata.var.reindex(new_var_names),
            X = scipy.sparse.csr_matrix((adata.shape[0], len(new_var_names))),
            obs = adata.obs,
            obsm=adata.obsm,
            uns=adata.uns
    )
    old_var_names = [x for x in new_var_names if x in adata.var_names]
    tmp_ad.X[:, new_var_names.isin(adata.var_names)] = adata[:, old_var_names].X
    return tmp_ad


# %%
query_hvg = reindex_adata(query, reference.var_names)

# %%
query_hvg

# %%
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
query.obs["predictions"] = vae_q.predict()


# %%
sc.pp.neighbors(query, use_rep="X_scANVI")
sc.tl.leiden(query)
sc.tl.umap(query)

# %%
sc.pl.umap(query, color=["predictions"])

# %%
adata_full = anndata.concat(reference, reindex_adata(query, reference.var_names))

# %%

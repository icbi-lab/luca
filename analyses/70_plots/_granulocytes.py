# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

# %%
import scanpy as sc
import pandas as pd

# %%
adata = sc.read_h5ad(
    "../../data/20_integrate_scrnaseq_data/28_annotate_cell_types_coarse_umap/adata_granulocytes.umap_leiden.h5ad"
)

# %%
pd.set_option("display.max_rows", 1000)

# %%
adata.obs.groupby(["dataset", "patient", "sample", "origin"], observed=True).size().reset_index(
    name="n_cells"
).sort_values("n_cells", ascending=False).to_csv("./granulocyte_count_per_patient.tsv", sep="\t")

# %%
sc.pl.umap(adata, color=["cell_type", "dataset"])

# %%
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)

# %%
sc.tl.pca(adata, use_highly_variable=True)

# %%
sc.pp.neighbors(adata)

# %%
sc.tl.umap(adata)

# %%
sc.pl.umap(adata, color=["cell_type", "dataset"])

# %%
adata

# %%

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
import threadpoolctl
import scipy
import numpy as np
import matplotlib.pyplot as plt
import scanpy_helpers as sh

# %%
threadpoolctl.threadpool_limits(20)

# %%
threadpoolctl.threadpool_info()

# %%
adata_unintegrated = sc.read_h5ad("../../data/20_build_atlas/integrate_datasets/21_merge_all/artifacts/merged_all.h5ad")

# %%
adata = sc.read_h5ad("../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad")

# %%
adata.X

# %%
adata.shape[0] * adata.shape[1]

# %%
adata_unintegrated = adata_unintegrated[adata.obs_names, :].copy()

# %%
sc.pp.normalize_total(adata_unintegrated, target_sum=1000)
sc.pp.log1p(adata_unintegrated)
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=6000)

# %%
sc.tl.pca(adata_unintegrated, svd_solver="arpack")

# %%
sc.pp.neighbors(adata_unintegrated)

# %%
sc.tl.umap(adata_unintegrated)

# %%
adata_unintegrated.obs["cell_type_coarse"] = adata.obs["cell_type_coarse"]
adata_unintegrated.obs["dataset"] = adata.obs["dataset"]
adata_unintegrated.obs["platform"] = adata.obs["platform"]

# %%
with plt.rc_context({"figure.figsize": (7, 7), "figure.dpi": 300}):
    sc.pl.umap(adata, color="cell_type_coarse")

# %%
adata.obs["platform"].cat.categories

# %%
sh.colors.set_scale_anndata(adata_unintegrated, "platform")

# %%
for c in ["dataset", "cell_type_coarse", "platform"]:
    with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": 300}):
        sc.pl.umap(adata_unintegrated, color=c, size=3, frameon=False)

# %%

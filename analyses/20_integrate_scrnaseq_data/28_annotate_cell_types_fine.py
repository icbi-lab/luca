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
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
from scanpy_helpers.annotation import AnnotationHelper
import scvi
import warnings
import numpy as np
from nxfvars import nxfvars

warnings.filterwarnings("ignore", category=FutureWarning)

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
input_dir = nxfvars.get("input_dir", "../../data/20_integrate_scrnaseq_data/27_annotate_cell_types_coarse_umap/")
main_adata = nxfvars.get("main_adata", "/home/sturm/Downloads/adata_cell_type_coarse.h5ad")
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(main_adata)

# %%
sc.pl.umap(adata, color="cell_type")

# %% [markdown]
# ## T cell subclustering

# %%
adata_t = sc.read_h5ad(f"{input_dir}/adata_t_cell.umap_leiden.h5ad")

# %%
ah.plot_umap(adata_t, filter_cell_type=["T cell", "NK"], cmap="inferno", size=5)

# %%
sc.pl.umap(adata_t, color="dataset")

# %%
adata_t.obs

# %%
sc.pl.umap(adata_t, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
cell_type_map = {
    "NK cell": []
}

# %% [markdown]
# ## Stromal cell subclustering

# %%
adata_stromal = sc.read_h5ad(f"{input_dir}/adata_t_cell.umap_leiden.h5ad")

# %%
adata_stromal.obs

# %%

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
import infercnvpy as cnv

warnings.filterwarnings("ignore", category=FutureWarning)

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
input_dir = nxfvars.get(
    "input_dir",
    "../../data/20_integrate_scrnaseq_data/28_annotate_cell_types_coarse_umap/",
)
main_adata = nxfvars.get(
    "main_adata",
    "../../data/20_integrate_scrnaseq_data/27_annotate_cell_types_coarse/artifacts/adata_cell_type_coarse.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(main_adata)

# %%
sc.pl.umap(adata, color="cell_type")

# %%
adata.obs["cell_type_coarse"] = adata.obs["cell_type"]

# %% [markdown]
# ## T cell subclustering

# %%
adata_t = sc.read_h5ad(f"{input_dir}/adata_t_cell.umap_leiden.h5ad")
adata_t.obs["leiden"] = adata_t.obs["leiden_1.00"]

# %%
ah.plot_umap(adata_t, filter_cell_type=["T cell", "NK", "Div"], cmap="inferno", size=1)

# %%
ah.plot_dotplot(adata_t, groupby="leiden_1.00")

# %%
sc.pl.umap(adata_t, color="dataset")

# %%
sc.pl.umap(adata_t, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
cell_type_map = {
    "B cell dividing": [17],
    "NK cell": [0, 18],
    "T cell dividing": [14],
    "T cell CD4": [2, 3, 5],
    "T cell CD8": [13, 1, 12, 8, 9, 11, 16, 19, 15],
    "T cell regulatory": [4, 6, 7],
    "T cell other": [10],
}

# %%
ah.annotate_cell_types(adata_t, cell_type_map)

# %%
ah.integrate_back(adata, adata_t)

# %% [markdown]
# ## Stromal cell subclustering

# %%
adata_stromal = sc.read_h5ad(f"{input_dir}/adata_stromal.umap_leiden.h5ad")

# %%
ah.plot_umap(
    adata_stromal,
    filter_cell_type=["Fibro", "muscle", "Peri", "Meso"],
    cmap="inferno",
    size=5,
)

# %%
adata_stromal.obs["leiden"] = adata_stromal.obs["leiden_0.50"]

# %%
sc.pl.umap(adata_stromal, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "Mesothelial": [8],
    "Pericyte": [3],
    "Smooth muscle cell": [6],
    "Fibroblast adventitial": [1],
    "Fibroblast alevolar": [2],
    "Fibroblast": [0, 9, 4, 5, 7],
}

# %%
ah.annotate_cell_types(adata_stromal, ct_map)

# %%
ah.integrate_back(adata, adata_stromal)

# %% [markdown]
# ## Myeloid compartment

# %%
adata_m = sc.read_h5ad(f"{input_dir}/adata_myeloid.umap_leiden.h5ad")

# %%
ah.plot_umap(
    adata_m, filter_cell_type=["Macro", "Mono", "DC", "Div"], cmap="inferno", size=2
)

# %%
ah.plot_umap(adata_m, filter_cell_type=["Alevolar"], cmap="inferno", size=2)

# %%
adata_m.obs["leiden"] = adata_m.obs["leiden_0.75"]

# %%
sc.pl.umap(adata_m, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.plot_dotplot(adata_m, groupby="leiden")

# %%
ct_map = {
    "myeloid dividing": [11],
    "DC mature/cDC 1": [14],
    "Macrophage FABP4+": [5, 0, 15, 9, 6, 10, 7],
    "Macrophage": [12, 1],
    "Monocyte": [4, 3, 8],
    "cDC2": [2, 16],
    "potential myeloid/epithelial doublets": [13],
}

# %%
ah.annotate_cell_types(adata_m, ct_map)

# %%
ah.integrate_back(adata, adata_m)

# %% [markdown]
# ## Epithelial compartment
#
# separate notebook

# %% [markdown]
# ## Write out results

# %%
adata.write_h5ad(f"{artifact_dir}/adata_annotated_fine.h5ad")

# %%

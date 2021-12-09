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
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/31_cell_types_coarse/by_cell_type/",
)
main_adata = nxfvars.get(
    "main_adata",
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/31_cell_types_coarse/artifacts/adata_cell_type_coarse.h5ad",
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
adata_t = sc.read_h5ad(f"{input_dir}/adata_cell_type_coarse_t_cell.umap_leiden.h5ad")
adata_t.obs["leiden"] = adata_t.obs["leiden_1.00"]

# %%
ah.plot_umap(adata_t, filter_cell_type=["T cell", "NK", "Div"], cmap="inferno", size=1)

# %%
ah.plot_dotplot(adata_t, groupby="leiden")

# %%
sc.pl.umap(adata_t, color="dataset")

# %%
sc.pl.umap(adata_t, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
cell_type_map = {
    "NK cell": [5, 3],
    "T cell dividing": [15],
    "T cell CD4": [1, 4, 0, 17, 16, 7],
    "T cell CD8": [13, 10, 8, 6, 12, 14, 9],
    "T cell regulatory": [2],
    "B cell dividing": [18],
    "other (T assoc.)": [11]
}

# %%
ah.annotate_cell_types(adata_t, cell_type_map)

# %%
ah.integrate_back(adata, adata_t)

# %% [markdown]
# ## Stromal cell subclustering

# %%
adata_stromal = sc.read_h5ad(
    f"{input_dir}/adata_cell_type_coarse_stromal.umap_leiden.h5ad"
)

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
    "Pericyte": [4],
    "Smooth muscle cell": [6],
    "Fibroblast adventitial": [0],
    "Fibroblast alevolar": [2],
    "Fibroblast": [1, 7, 5, 3, 9],
}

# %%
ah.annotate_cell_types(adata_stromal, ct_map)

# %%
ah.integrate_back(adata, adata_stromal)

# %% [markdown]
# ## Myeloid compartment

# %%
adata_m = sc.read_h5ad(f"{input_dir}/adata_cell_type_coarse_myeloid.umap_leiden.h5ad")

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
    "myeloid dividing": [16, 11],
    "Macrophage FABP4+": [0, 2, 9, 7],
    "Macrophage": [3, 6, 15, 12, 5, 10],
    "Monocyte": [1, 14, 8, 17],
    "cDC2": [4],
    "other (myeloid-assoc.)": [18],
    "cDC1/mature": [13]
}

# %%
ah.annotate_cell_types(adata_m, ct_map)

# %%
adata_dc = adata_m[adata_m.obs["leiden"] == "13", :]

# %%
ah.reprocess_adata_subset_scvi(adata_dc, leiden_res=0.5)

# %%
ah.plot_umap(adata_dc, filter_cell_type=["DC"], cmap="inferno")

# %%
sc.pl.umap(adata_dc, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(adata_dc, {"DC mature": [2, 0, 8, 9], "cDC1": [5, 3, 6, 4, 1, 7]})

# %%
ah.integrate_back(adata_m, adata_dc)

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

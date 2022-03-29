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
import warnings
import numpy as np
from nxfvars import nxfvars
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)

# %%
sc.set_figure_params(figsize=(3, 3))

# %%
# based on curated marker genes from several studies
ah = AnnotationHelper()

# %%
# based on Human Lung Cell Atlas
ah2 = AnnotationHelper(
    markers=pd.read_csv(
        "https://docs.google.com/spreadsheets/d/1beW-9oeM31P50NFvNLVvsdXlfh_tOjsmIrnwV2ZlxDU/gviz/tq?tqx=out:csv&sheet=lung_hlca"
    )
)

# %%
input_dir = nxfvars.get(
    "input_dir",
    "../../data/20_build_atlas//annotate_datasets/31_cell_types_coarse/by_cell_type/",
)
main_adata = nxfvars.get(
    "main_adata",
    "../../data/20_build_atlas/annotate_datasets/31_cell_types_coarse/artifacts/adata_cell_type_coarse.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(main_adata)

# %%
sc.pl.umap(adata, color="cell_type")

# %%
adata.obs["cell_type_coarse"] = adata.obs["cell_type"]

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# ## B cells

# %%
adata_b = sc.read_h5ad(f"{input_dir}/adata_cell_type_coarse_b_cell.umap_leiden.h5ad")

# %%
ah2.score_cell_types(adata_b)

# %%
ah.plot_umap(
    adata_b,
    filter_cell_type=["b cell", "div"],
    cmap="inferno",
    size=5,
)

# %%
adata_b.obs["leiden"] = adata_b.obs["leiden_0.50"]

# %%
ah.plot_dotplot(adata_b)

# %%
ah2.plot_dotplot_scores(adata_b)

# %%
sc.pl.umap(adata_b, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(adata_b, cell_type_map={"B cell": [0, 1, 2, 3, 4, 5, 6], "potential doublets": [7]})

# %%
ah.integrate_back(adata, adata_b)

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true
# ## Endothelial cell subclustering

# %%
adata_endo = sc.read_h5ad(
    f"{input_dir}/adata_cell_type_coarse_endothelial_cell.umap_leiden.h5ad"
)
adata_endo.obs["leiden"] = adata_endo.obs["leiden_0.75"]

# %%
sc.pl.umap(adata_endo, color=["origin", "condition", "dataset"])

# %%
ah.plot_umap(adata_endo, filter_cell_type=["Endo", "Div"], cmap="inferno", size=2)

# %%
ah2.score_cell_types(adata_endo)

# %%
ah2.plot_umap_scores(adata_endo, filter_cell_type=["EC", "Endo"])

# %%
ah2.plot_dotplot_scores(adata_endo)

# %%
with plt.rc_context({"figure.figsize": (5, 5)}):
    sc.pl.umap(adata_endo, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(
    adata_endo,
    {
        "Endothelial cell lymphatic": [2],
        "Endothelial cell venous": [0, 3, 4, 10],
        "Endothelial cell arterial": [7],
        "Endothelial cell capillary": [9, 8, 6, 1],
        "Endothelial cell": [5],
        "potential doublets": [11],
    },
)

# %%
ah2.plot_dotplot_scores(
    adata_endo, filter_cell_type=["EC", "B cell"], groupby="cell_type"
)

# %%
ah.integrate_back(adata, adata_endo)

# %% [markdown]
# ## Epithelial cells
# --> separate notebook

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# ## Granulocytes

# %%
adata_g = sc.read_h5ad(f"{input_dir}/adata_cell_type_coarse_granulocytes.umap_leiden.h5ad")
adata_g.obs["leiden"] = adata_g.obs["leiden_1.00"]

# %%
ah.plot_umap(adata_g, cmap="inferno", filter_cell_type=["Neutro", "MDSC", "Granul", "div"])

# %%
ah.plot_dotplot(adata_g)

# %%
ah2.score_cell_types(adata_g)

# %%
ah2.plot_dotplot_scores(adata_g)

# %%
sc.pl.umap(adata_g, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(adata_g, cell_type_map={
    "Neutrophils": [2, 0, 10, 8, 4, 12, 9, 3, 1, 7, 6, 5, 11],
    "potential doublets": [13]
})

# %%
ah.integrate_back(adata, adata_g)

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true
# ## Mast cells

# %%
adata_mast = sc.read_h5ad(
    f"{input_dir}/adata_cell_type_coarse_mast_cell.umap_leiden.h5ad"
)

# %%
ah2.score_cell_types(adata_mast)

# %%
ah.plot_umap(
    adata_mast,
    filter_cell_type=["mast", "div"],
    cmap="inferno",
    size=5,
)

# %%
sc.pl.umap(adata_mast, color=["dataset"])

# %%
adata_mast.obs["leiden"] = adata_mast.obs["leiden_1.00"]

# %%
ah.plot_dotplot(adata_mast)

# %%
ah2.plot_dotplot_scores(adata_mast)

# %%
sc.pl.umap(adata_mast, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(
    adata_mast,
    cell_type_map={
        "Mast cell": [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 16],
        "potential doublets": [9, 17, 15]
    },
)

# %% tags=[]
ah.integrate_back(adata, adata_mast)

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true
# ## Myeloid compartment

# %%
adata_m = sc.read_h5ad(f"{input_dir}/adata_cell_type_coarse_myeloid.umap_leiden.h5ad")

# %%
ah.plot_umap(
    adata_m, filter_cell_type=["Macro", "Mono", "DC", "Div"], cmap="inferno", size=2
)

# %%
adata_m.obs["log_counts"] = np.log1p(adata_m.obs["total_counts"])

# %%
sc.pl.umap(adata_m, color=["log_counts", "n_genes_by_counts"])

# %%
ah2.score_cell_types(adata_m)

# %%
ah2.plot_umap_scores(adata_m, filter_cell_type=["DC", "macro", "Mono"])

# %%
adata_m.obs["leiden"] = adata_m.obs["leiden_0.75"]

# %%
ah2.plot_dotplot_scores(adata_m)

# %%
sc.pl.umap(adata_m, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.plot_dotplot(adata_m, groupby="leiden")

# %%
ct_map = {
    "myeloid dividing": [11],
    "Macrophage alveolar": [0, 2, 9, 7],
    "Macrophage": [3, 6, 15, 12, 5, 10],
    "Monocyte classical": [1, 17, 8],
    "Monocyte non-classical": [14],
    "cDC2": [4],
    "cDC1/mature": [13],
    "potential epithelial/myeloid doublets": [16],
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
ah.annotate_cell_types(
    adata_dc, {"DC mature": [2, 0, 8, 9], "cDC1": [5, 3, 6, 4, 1, 7]}
)

# %%
ah.integrate_back(adata_m, adata_dc)

# %%
ah.integrate_back(adata, adata_m)

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true
# ## Plasma cells

# %%
adata_p = sc.read_h5ad(
    f"{input_dir}/adata_cell_type_coarse_plasma_cell.umap_leiden.h5ad"
)

# %%
adata_p.obs["leiden"] = adata_p.obs["leiden_0.50"]

# %%
ah2.score_cell_types(adata_p)

# %%
ah.plot_umap(adata_p, filter_cell_type=["B cell", "Plasma"], cmap="inferno", size=2)

# %%
ah.plot_dotplot(adata_p)

# %%
ah2.plot_dotplot_scores(adata_p)

# %%
adata_p.obs["log_counts"] = np.log1p(adata_p.obs["total_counts"])

# %%
sc.pl.umap(adata_p, color=["log_counts", "n_genes_by_counts"])

# %%
sc.pl.umap(adata_p, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(
    adata_p,
    cell_type_map={
        "Plasma cell": [6, 7, 2, 5, 0, 3, 1],
        "Plasma cell dividing": [8],
        "unknown": [9, 10],
        "potential plasma/myeloid doublets": [4],
    },
)

# %%
ah.integrate_back(adata, adata_p)

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# ## Stromal cell subclustering

# %%
adata_stromal = sc.read_h5ad(
    f"{input_dir}/adata_cell_type_coarse_stromal.umap_leiden.h5ad"
)

# %%
ah.plot_umap(
    adata_stromal,
    filter_cell_type=["Fibro", "muscle", "Peri", "Meso", "div"],
    cmap="inferno",
    size=5,
)

# %%
ah2.score_cell_types(adata_stromal)

# %%
ah2.plot_umap_scores(
    adata_stromal, filter_cell_type=["fibro", "muscle", "SM", "Peri", "Fibro", "Meso"]
)

# %%
adata_stromal.obs["leiden"] = adata_stromal.obs["leiden_0.50"]

# %%
ah2.plot_dotplot_scores(
    adata_stromal
)

# %%
sc.pl.umap(adata_stromal, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "Mesothelial": [8],
    "Pericyte": [4],
    "Smooth muscle cell": [6],
    "Fibroblast adventitial": [0],
    "Fibroblast alveolar": [2],
    "Fibroblast peribronchial": [1],
    "stromal": [5, 3, 7],
    "stromal dividing": [9],
}

# %%
ah.annotate_cell_types(adata_stromal, ct_map)

# %%
ah.integrate_back(adata, adata_stromal)

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true
# ## T cell subclustering

# %%
adata_t = sc.read_h5ad(f"{input_dir}/adata_cell_type_coarse_t_cell.umap_leiden.h5ad")
adata_t.obs["leiden"] = adata_t.obs["leiden_1.00"]

# %%
ah.plot_umap(adata_t, filter_cell_type=["T cell", "NK", "Div"], cmap="inferno", size=1)

# %%
ah.plot_dotplot(adata_t, groupby="leiden")

# %% [markdown]
# ### What's cluster 11? 

# %%
sc.tl.rank_genes_groups(adata_t, "leiden", groups=["11"], reference="rest")

# %%
sc.pl.rank_genes_groups_dotplot(adata_t, dendrogram=False)

# %%
adata_t.obs["log_counts"] = np.log1p(adata_t.obs["total_counts"])

# %%
sc.pl.umap(adata_t, color=["log_counts", "n_genes_by_counts"])

# %%
sc.pl.umap(adata_t, color="dataset")

# %% [markdown]
# Cluster 11 comes from a single dataset, has very few genes and few counts. It does not express any markers known to us. I suspect these are empty droplets that passed through QC. 

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
    "potentially empty droplets": [11],
}

# %%
ah.annotate_cell_types(adata_t, cell_type_map)

# %%
ah.integrate_back(adata, adata_t)

# %% [markdown]
# ## pDC

# %%
adata_pdc = sc.read_h5ad(f"{input_dir}/adata_cell_type_coarse_pdc.umap_leiden.h5ad")
adata_pdc.obs["leiden"] = adata_pdc.obs["leiden_1.00"] 

# %%
ah.plot_umap(adata_pdc, filter_cell_type=["div", "pdc"], cmap="inferno")

# %%
ah.plot_dotplot(adata_pdc)

# %%
ah2.score_cell_types(adata_pdc)

# %%
ah2.plot_dotplot_scores(adata_pdc)

# %%
sc.pl.umap(adata_pdc, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(adata_pdc, cell_type_map={"pDC": range(11), "potential doublets": [11]})

# %%
ah.integrate_back(adata, adata_pdc)

# %% [markdown]
# ## Write out results

# %%
adata.write_h5ad(f"{artifact_dir}/adata_annotated_fine.h5ad")

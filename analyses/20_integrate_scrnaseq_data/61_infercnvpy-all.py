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
import infercnvpy as cnv
import scanpy as sc
from scanpy_helpers.annotation import AnnotationHelper
import scvi
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

# %%
ah = AnnotationHelper()

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
adata = sc.read_h5ad(
    "../../data/50_integrate_scrnaseq_data/52_run_scanvi/integrated_merged_all_all_genes.h5ad"
)

# %%
adata.shape

# %%
scvi_model = scvi.model.SCVI.load(
    "../../data/50_integrate_scrnaseq_data/52_run_scanvi/scvi_model_merged_all_all_genes//",
    adata=adata,
)

# %%
adata.obs["batch"][0]

# %%
# %%time
solo = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch="Adams_Kaminski_2020_COPD_Adams_Kaminski_2020_COPD_001C")

# %%
adata_infercnv = sc.AnnData(
    X=adata.raw.X,
    var=adata.raw.var,
    obs=adata.obs,
    obsm=adata.obsm,
    uns=adata.uns,
    obsp=adata.obsp,
    layers={"imputed": scvi_model.get_normalized_expression(library_size=10000)}
)

# %%
cnv.io.genomic_position_from_gtf(
    "/data/genomes/hg38/annotation/gencode/gencode.v33.primary_assembly.annotation.gtf",
    adata_infercnv,
)

# %%
sc.pl.umap(adata, color=["dataset"])

# %%
sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=30)
sc.tl.umap(adata)

# %%
sc.pl.umap(adata, color=["dataset"])

# %%
# %%time
cnv.tl.infercnv(
    adata_infercnv,
    reference_key="cell_type_predicted",
    reference_cat=[
        "T cell CD4",
        "T cell CD8",
        "Monocyte conventional",
        "Macrophage",
        "Mast cell",
        "Plasma cell",
        "B cell",
        "pDC",
    ],
    window_size=250,
)

# %%
adata = adata_infercnv

# %%
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.umap(adata)

# %%
cnv.tl.leiden(adata, resolution=1.5)

# %%
cnv.tl.cnv_score(adata)

# %%
cnv.pl.umap(adata, color=["cell_type", "cnv_leiden", "cnv_score"])

# %%
sc.pl.umap(adata, color="cell_type_predicted")

# %%
sc.tl.leiden(adata)

# %%
adata.write_h5ad("/home/sturm/Downloads/tmp_adata_61.h5ad")

# %%
adata = sc.read_h5ad("/home/sturm/Downloads/tmp_adata_61.h5ad")

# %%
sc.pl.umap(adata, color="leiden")

# %%
ah.plot_umap(adata)

# %%
ah.plot_dotplot(adata)

# %%
sc.pl.umap(adata, color=["cnv_score", "cnv_leiden", "origin"])

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "B cell": [14],
    "Ciliated": [18],
    "Endothelial cell": [6],
    "Endothelial cell lymphatic": [22],
    "Epithelial cell": [15, 24, 21, 19, 3, 37, 13, 31, 20, 12, 32],
    "Mast cell": [21],
    "Myeloid": [4, 0, 10, 1, 9, 29, 34, 28, 44, 26, 5, 8, 25, 40, 27, 38],
    "Stromal": [11],
    "Plasma cell": [23],
    "T cell": [2, 36, 17, 16, 35, 4, 7],
    "pDC": [30],
}

# %%
ah.annotate_cell_types(adata, ct_map)

# %%
sc.pl.umap(adata, color="cell_type", groups=["other"])

# %%
cnv.pl.umap(adata, color="cell_type")

# %%
adata_t = adata[adata.obs["cell_type"] == "T cell", :].copy()
ah.reprocess_adata_subset_scvi(adata_t, use_rep="X_scANVI", leiden_res=1.5)

# %%
ah.plot_umap(adata_t, filter_cell_type=["T cell", "NK", "Div"])

# %%
sc.pl.umap(adata_t, color=["KLRB1"])

# %%
sc.pl.umap(adata_t, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.plot_dotplot(adata_t)

# %%
ct_mat = {
    "T cell dividing": [16],
    "NK cell": [1, 0, 17],
    "T cell CD8": [6, 7, 22, 14, 10, 23, 3, 13, 24, 25],
    "T reg": [5],
    "T cell CD4": [9, 2, 19, 4, 8, 15],
    "other T assoc.": [11, 12, 21, 18, 20, 27, 26],
}

# %%
ah.annotate_cell_types(adata_t, ct_mat)

# %%
ah.integrate_back(adata, adata_t)

# %%
adata_m = adata[adata.obs["cell_type"] == "Myeloid", :].copy()
ah.reprocess_adata_subset_scvi(adata_m, use_rep="X_scANVI", n_neighbors=15)

# %%
sc.pl.umap(adata_m, color=["CD1A", "CD207", "MARCO", "CD163", "SLAMF9"])

# %%
ah.plot_umap(adata_m, filter_cell_type=["Macro", "Mono", "DC", "Div", "MDSC"])

# %%
ah.plot_dotplot(adata_m)

# %%
sc.pl.umap(adata_m, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "mDC mature": [17],
    "cDC2": [5, 23],
    "myeloid dividing": [15, 21],
    "Macrophage FABP4+": [19, 3, 0, 2, 24, 10],
    "Macrophage": [8, 4, 1, 13, 18],
    "Monocyte conventional": [6],
    "Monocyte non-conventional": [11],
    "cDC1": [16],
    "MDSC": [12, 20, 22],
    "myeloid other": [9, 7, 14]
}

# %%
ah.annotate_cell_types(adata_m, ct_map)

# %%
ah.integrate_back(adata, adata_m)

# %%
adata_s = adata[adata.obs["cell_type"] == "Stromal", :].copy()
ah.reprocess_adata_subset_scvi(
    adata_s, leiden_res=1, use_rep="X_scANVI", n_neighbors=15
)

# %%
ah.plot_umap(
    adata_s, filter_cell_type=["Fibro", "muscle", "Peri", "Meso", "Endo", "Div"]
)

# %%
sc.pl.umap(adata_s, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "Fibroblast": [2, 16, 7, 8, 17, 10, 15, 13, 19],
    "Fibroblast adventitial": [6, 0, 9],
    "Fibroblast alevolar": [4, 18, 1],
    "Smooth muscle cell": [3, 12],
    "Mesothelial": [11],
    "Pericyte": [5],
    "stromal other": [14],
}

# %%
ah.annotate_cell_types(adata_s, ct_map)

# %%
ah.integrate_back(adata, adata_s)

# %%
adata_epi = adata[adata.obs["cell_type"] == "Epithelial cell", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_epi, use_rep="X_scANVI")

# %%
ah.reprocess_adata_subset_cnv(adata_epi)

# %%
cnv.pl.umap(
    adata_epi,
    color=["cnv_leiden", "cnv_score", "EPCAM"],
    legend_loc="on data",
    legend_fontoutline=2,
)

# %%
cnv.tl.infercnv(adata_epi, reference_key="cnv_leiden", reference_cat=["4", "2", "10"])

# %%
cnv.pl.chromosome_heatmap(adata_epi)

# %%
sc.pl.umap(
    adata_epi, color=["cnv_leiden", "cnv_score", "EPCAM", "dataset", "origin"], ncols=2
)

# %%
ah.plot_umap(
    adata_epi,
    filter_cell_type=[
        "Alevolar",
        "Basal",
        "Club",
        "Dividing",
        "Goblet",
        "Ionocyte",
        "Mesothelial",
        "Suprabasal",
        "Epi",
    ],
)

# %%
ah.plot_dotplot(adata_epi)

# %%
sc.pl.umap(adata_epi, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ct_map = {
    "Epithelial cell (malignant)": [12, 25, 28, 14, 17, 13, 21, 32],
    "Alevolar cell type 1": [0, 26, 18],
    "Alevolar cell type 2": [16, 3, 4, 2, 9, 20, 27, 22, 1, 10, 8, 11, 15, 29, 30],
    "Club/Goblet": [7, 5, 23],
    "Epithelial cell other (benign)": [6, 24, 19, 31],
}

# %%
ah.annotate_cell_types(adata_epi, ct_map)

# %%
ah.integrate_back(adata, adata_epi)

# %%
adata.obs["is_malignant"] = [
    "malignant" if x else "non-malignant"
    for x in adata.obs["cell_type"].str.contains("mali")
]

# %%
sc.pl.umap(adata, color=["is_malignant", "cnv_score", "cnv_leiden", "origin", "dataset"])

# %%
adata.obs["condition_origin"] = [
    f"{condition}_{origin}"
    for condition, origin in zip(adata.obs["condition"], adata.obs["origin"])
]

# %%
adata.write_h5ad("../../data/60_infercnv/all-annotated-integrated.h5ad")

# %%
# !mkdir -p ../../tables/single_cell_annotations/all_integrated

# %%
adata_scvi = sc.AnnData(
    X=adata.obsm["X_scANVI"], obs=adata.obs.loc[:, ["cell_type", "cnv_leiden"]]
)
adata_scvi.write_h5ad(
    "../../tables/single_cell_annotations/all_integrated/adata_scvi.h5ad",
    compression="gzip",
    compression_opts=9,
)

# %%
adata.obs.loc[:, ["cell_type", "cnv_leiden"]].to_csv(
    "../../tables/single_cell_annotations/all_integrated/cell_type_annotations.csv"
)

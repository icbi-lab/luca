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
import pandas as pd
import scanpy_helpers as sh
import scanpy as sc
import itertools

sc.settings.set_figure_params(figsize=(5, 5))

# %%
de_res = pd.read_csv(
    "../../data/30_downstream_analyses/de_analysis/tumor_cell_types/de_deseq2/core_atlas_tumor_cell_types_core_atlas_tumor_cells_DESeq2_result.tsv",
    sep="\t",
)

# %%
adata = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/33_cell_types_epi/artifacts/adata_tumor.h5ad"
)

# %%
de_res

# %%
sc.pl.umap(adata, color=["cell_type"])

# %%
adata = adata[adata.obs["cell_type"].str.startswith("Tumor"), :].copy()

# %%
top_markers = (
    de_res.loc[lambda x: (x["log2FoldChange"] > 0) & (x["baseMean"] > 20)]
    .groupby("comparison")
    .apply(lambda x: x.head(10))
    .reset_index(drop=True)
)

# %%
pd.set_option("display.max_rows", None)

# %%
de_res.loc[lambda x: x["log2FoldChange"] > 0].to_csv("/home/sturm/Downloads/tumor_clusters_all_genes.tsv", sep="\t")

# %%
top_markers.to_csv("/home/sturm/Downloads/tumor_clusters_top_markers.tsv", sep="\t")

# %%
top_markers

# %%
sh.compare_groups.pl.plot_lm_result_altair(
    de_res.loc[lambda x: x["gene_id"].isin(top_markers["gene_id"].values)],
    x="gene_id",
    y="comparison",
    p_col="padj",
    color="log2FoldChange",
    cluster=False,
    value_max=5,
    order=top_markers["gene_id"].tolist()
)

# %%
pb = sh.pseudobulk.pseudobulk(adata, groupby=["patient", "cell_type"])

# %%
sc.pp.normalize_total(pb, target_sum=1e6)
sc.pp.log1p(pb)

# %%
sc.pl.matrixplot(
    pb,
    groupby="cell_type",
    var_names=top_markers.groupby("comparison")
    .apply(lambda x: x["gene_id"].tolist())
    .to_dict(),
)

# %%
tumor_markers = {
    # "Epithelial": ["EPCAM", "B2M"],
    # "Alveolar cell type 1": ["AGER", "CLDN18"],
    # "Alveolar cell type 2": ["SFTPC", "SFTPB", "SFTPA1"],
    # "Ciliated": ["PIFO", "FOXJ1", "HYDIN", "CFAP299"],
    # "Club": ["SCGB3A1", "SCGB3A2"],
    "LUAD": ["CD24", "MUC1", "NAPSA", "NKX2-1", "KRT7", "MSLN"],
    "LUSC": ["KRT5", "KRT6A", "TP63", "NTRK2", "SOX2", "KRT17"],
    "NE": ["CHGA", "SYP", "NCAM1", "TUBA1A"],
    "EMT": ["VIM", "SERPINE1", "CDH1"],
    "mitotic": ["MKI67", "TOP2A"],
    "undifferentiated": ["TACSTD2", "AGR2"],
}

# %%
known_markers = list(itertools.chain.from_iterable(tumor_markers.values()))

# %%
sh.compare_groups.pl.plot_lm_result_altair(
    de_res.loc[lambda x: x["gene_id"].isin(known_markers)],
    x="gene_id",
    y="comparison",
    p_col="padj",
    color="log2FoldChange",
    cluster=False,
    value_max=5,
    order=known_markers,
)

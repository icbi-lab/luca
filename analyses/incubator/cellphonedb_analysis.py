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
from nxfvars import nxfvars
import altair as alt
from toolz.functoolz import pipe
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy_helpers as sh
import statsmodels.formula.api as smf
from tqdm.contrib.concurrent import process_map
import itertools
import progeny
import dorothea
import scipy.stats
from threadpoolctl import threadpool_limits
import warnings

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
adata_primary_tumor = adata[(adata.obs["origin"] == "tumor_primary")].copy()

# %%
adata_n = sc.read_h5ad(
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/adata_neutrophil_clusters.h5ad"
    # "/home/sturm/Downloads/adata_neutrophil_clusters.h5ad"
)

# %%
immune_cells = [
    "B cell",
    "cDC1",
    "cDC2",
    "DC mature",
    "Macrophage",
    "Macrophage alveolar",
    "Mast cell",
    "Monocyte",
    "Neutrophils",
    "NK cell",
    "pDC",
    "Plasma cell",
    "T cell CD4",
    "T cell CD8",
    "T cell regulatory",
]

# %% [markdown]
# ## Neutrophils

# %%
pb_adata = sh.pseudobulk.pseudobulk(
    adata_n, groupby=["dataset", "patient", "cell_type"]
)

# %%
sc.pp.normalize_total(pb_adata, target_sum=1e6)
sc.pp.log1p(pb_adata)

# %%
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    de_genes = sh.compare_groups.lm.test_lm(
        pb_adata, "~ C(cell_type, Sum) + patient", groupby="cell_type", n_jobs=16
    )

# %%
# de_res_tumor_cells = (
#     pd.concat(
#         [
#             pd.read_csv(
#                 f"../../data/30_downstream_analyses/de_analysis/{k}_desert/de_deseq2/adata_primary_tumor_tumor_cells_DESeq2_result.tsv",
#                 sep="\t",
#             ).assign(group=k.upper())
#             for k in "tmb"
#         ]
#     )
#     .fillna(1)
#     .pipe(sh.util.fdr_correction)
#     .drop("gene_id", axis="columns")
#     .rename(columns={"gene_id.1": "gene_id"})
# )

# %% [markdown]
# ## LUAD vs LUSC

# %%
cpdb = pd.read_csv("../../tables/cellphonedb_2022-04-06.tsv", sep="\t")

# %%
cpdba = sh.cell2cell.CpdbAnalysis(
    cpdb,
    adata_primary_tumor,
    pseudobulk_group_by=["patient"],
    cell_type_column="cell_type_major",
)

# %%
de_res_tumor_cells = (
    pd.read_csv(
        "../../data/30_downstream_analyses/de_analysis/luad_lusc/de_deseq2/adata_primary_tumor_tumor_cells_DESeq2_result.tsv",
        sep="\t",
    )
    .fillna(1)
    .pipe(sh.util.fdr_correction)
    .drop("gene_id", axis="columns")
    .rename(columns={"gene_id.1": "gene_id"})
    .assign(group="LUSC")
)

# %%
cpdb_res = cpdba.significant_interactions(de_res_tumor_cells, max_pvalue=0.1)

# %%
cpdb_res

# %%
top_genes = (
    cpdb_res.loc[:, ["source_genesymbol", "fdr"]]
    .drop_duplicates()
    .sort_values("fdr")["source_genesymbol"][:30]
    .tolist()
)

# %%
cpdba.plot_result(
    cpdb_res.loc[
        lambda x: x["cell_type_major"].isin(immune_cells)
        & x["source_genesymbol"].isin(top_genes)
    ],
    title="LUAD vs LUSC: tumor cells",
    aggregate=False,
    cluster="heatmap",
)

# %% [markdown]
# ---

# %%
expressed_genes.to_csv("/home/sturm/Downloads/expressed_genes.csv")
de_res_tumor_cells.loc[
    lambda x: x["gene_id"].isin(significant_interactions["source_genesymbol"])
].assign(log2FoldChange=lambda x: np.clip(x["log2FoldChange"], -5, 5)).to_csv(
    "/home/sturm/Downloads/de_tumor_ligands.csv"
)

# %%

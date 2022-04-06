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
cpdb = pd.read_csv("../../tables/cellphonedb_2022-04-06.tsv", sep="\t")

# %%
cpdb

# %%
significant_genes = de_res_tumor_cells.loc[
    lambda x: x["fdr"] < 0.001, "gene_id"
].unique()

# %%
significant_interactions = cpdb.loc[
    lambda x: x["source_genesymbol"].isin(significant_genes)
]

# %%
pb_fracs = sh.pseudobulk.pseudobulk(
    adata_primary_tumor,
    groupby=["dataset", "patient", "cell_type_major"],
    aggr_fun=lambda x, axis: np.sum(x > 0, axis) / x.shape[axis],
)

# %%
fractions_expressed = sh.pseudobulk.pseudobulk(
    pb_fracs, groupby="cell_type_major", aggr_fun=np.mean
)
fractions_expressed.obs.set_index("cell_type_major", inplace=True)

# %%
expressed_genes = (
    fractions_expressed.to_df()
    .melt(ignore_index=False, value_name="fraction_expressed")
    .loc[lambda x: x["fraction_expressed"] >= 0.1]
    .reset_index()
)

# %%
heatmap_df = (
    expressed_genes.merge(
        significant_interactions.loc[:, ["source_genesymbol", "target_genesymbol"]],
        left_on="variable",
        right_on="target_genesymbol",
    )
    .drop(columns=["variable"])
    .groupby(["cell_type_major", "source_genesymbol"])
    .agg(
        n=("target_genesymbol", len),
        max_frac_expressed=("fraction_expressed", np.max),
    )
).reset_index()

# %%
p1 = (
    de_res_tumor_cells.loc[
        lambda x: x["gene_id"].isin(significant_interactions["source_genesymbol"])
        & x["gene_id"].isin(heatmap_df["source_genesymbol"])
    ]
    .assign(log2FoldChange=lambda x: np.clip(x["log2FoldChange"], -5, 5))
    .pipe(
        sh.compare_groups.pl.plot_lm_result_altair,
        title="Differential genes (tumor cells)",
        color="log2FoldChange",
        x="gene_id",
        configure=lambda x: x,
    )
)
p1

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

# %%
(
    p1
    & alt.Chart(heatmap_df.loc[lambda x: x["cell_type_major"].isin(immune_cells)])
    .mark_circle()
    .encode(
        x=alt.X(
            "source_genesymbol",
            axis=alt.Axis(grid=True),
            # scale=alt.Scale(
            #     domain=np.unique(significant_interactions["source_genesymbol"].values)
            # ),
        ),
        y=alt.Y("cell_type_major", axis=alt.Axis(grid=True)),
        size=alt.Size("n:N", scale=alt.Scale(domain=[1, 2], range=[80, 140])),
        color=alt.Color("max_frac_expressed", scale=alt.Scale(scheme="cividis")),
    )
).resolve_scale(size="independent", color="independent", x="shared").configure_mark(
    opacity=1
)

# %%
heatmap_df

# %%
expressed_genes.to_csv("/home/sturm/Downloads/expressed_genes.csv")

# %%
de_res_tumor_cells.loc[
    lambda x: x["gene_id"].isin(significant_interactions["source_genesymbol"])
].assign(log2FoldChange=lambda x: np.clip(x["log2FoldChange"], -5, 5)).to_csv(
    "/home/sturm/Downloads/de_tumor_ligands.csv"
)

# %%

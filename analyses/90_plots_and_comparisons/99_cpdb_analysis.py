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
    "main_adata",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
path_adata_n = nxfvars.get(
    "adata_n",
    "../../data/30_downstream_analyses/neutrophils/subclustering//artifacts/adata_neutrophil_clusters.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")
path_cpdb = nxfvars.get("path_cpdb", "../../tables/cellphonedb_2022-04-06.tsv")
deseq2_path_prefix = nxfvars.get(
    "deseq2_path_prefix",
    "../../data/30_downstream_analyses/de_analysis/{comparison}/de_deseq2",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
adata_n = sc.read_h5ad(path_adata_n)

# %%
cpdb = pd.read_csv(path_cpdb, sep="\t")

# %%
adata_primary_tumor = adata[(adata.obs["origin"] == "tumor_primary")].copy()

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
cpdb

# %%
cpdba = sh.cell2cell.CpdbAnalysis(
    cpdb,
    adata_primary_tumor,
    pseudobulk_group_by=["patient"],
    cell_type_column="cell_type_major",
)

# %%
cpdba_luad = sh.cell2cell.CpdbAnalysis(
    cpdb,
    adata_primary_tumor[adata_primary_tumor.obs["condition"] == "LUAD", :].copy(),
    pseudobulk_group_by=["patient"],
    cell_type_column="cell_type_major",
)

# %%
cpdba_lusc = sh.cell2cell.CpdbAnalysis(
    cpdb,
    adata_primary_tumor[adata_primary_tumor.obs["condition"] == "LUSC", :].copy(),
    pseudobulk_group_by=["patient"],
    cell_type_column="cell_type_major",
)

# %%
cpdba_coarse = sh.cell2cell.CpdbAnalysis(
    cpdb,
    adata_primary_tumor,
    pseudobulk_group_by=["patient"],
    cell_type_column="cell_type_coarse",
)

# %% [markdown] tags=[]
# # LUAD vs LUSC

# %%
de_res_tumor_cells_luad_lusc = (
    pd.read_csv(
        (
            deseq2_path_prefix
            + "/{comparison}_primary_tumor_adata_primary_tumor_tumor_cells_DESeq2_result.tsv"
        ).format(comparison="luad_lusc"),
        sep="\t",
    )
    .fillna(1)
    .pipe(sh.util.fdr_correction)
    .drop("gene_id", axis="columns")
    .rename(columns={"gene_id.1": "gene_id"})
    .assign(group="LUAD")
)

# %%
cpdb_res = cpdba_coarse.significant_interactions(
    de_res_tumor_cells_luad_lusc, max_pvalue=0.1
)
cpdb_res.to_csv(f"{artifact_dir}/cpdb_luad_lusc_coarse.csv")

# %% [markdown]
# ### LUSC

# %%
sc.pl.dotplot(
    adata_primary_tumor[adata_primary_tumor.obs["cell_type_major"] == "Tumor cells", :],
    var_names=["PGF"],
    groupby="condition",
)

# %%
cpdb_res = cpdba_lusc.significant_interactions(
    de_res_tumor_cells_luad_lusc, max_pvalue=0.1
)
top_genes = (
    cpdb_res.loc[:, ["source_genesymbol", "fdr"]]
    .drop_duplicates()
    .sort_values("fdr")["source_genesymbol"][:30]
    .tolist()
)
# cpdb_res.to_csv(f"{artifact_dir}/cpdb_luad_lusc.csv")
cpdba_lusc.plot_result(
    cpdb_res.loc[
        lambda x: x["cell_type_major"].isin(immune_cells)
        & x["source_genesymbol"].isin(top_genes)
    ],
    title="LUAD vs LUSC: tumor cells, top 30 DE ligands",
    aggregate=False,
    cluster="heatmap",
    label_limit=80,
)

# %% [markdown]
# # Patient stratification

# %%
de_res_tumor_cells_patient_strat = (
    pd.concat(
        [
            pd.read_csv(
                (
                    deseq2_path_prefix
                    + "/{k}_desert_primary_tumor_adata_primary_tumor_tumor_cells_DESeq2_result.tsv"
                ).format(k=k),
                sep="\t",
            ).assign(group=k.upper())
            for k in "tmb"
        ]
    )
    .fillna(1)
    .pipe(sh.util.fdr_correction)
    .drop("gene_id", axis="columns")
    .rename(columns={"gene_id.1": "gene_id"})
)

# %%
cpdb_res = cpdba.significant_interactions(
    de_res_tumor_cells_patient_strat, max_pvalue=0.1
)
cpdb_res.to_csv(f"{artifact_dir}/cpdb_patient_stratification.csv")
cpdba.plot_result(
    cpdb_res.loc[lambda x: x["cell_type_major"].isin(immune_cells)],
    title="Patient stratification: tumor cells, ligand FDR < 0.1",
    aggregate=False,
    cluster="heatmap",
)

# %% [markdown]
# # Neutrophils

# %%
pb_adata = sh.pseudobulk.pseudobulk(
    adata_n, groupby=["dataset", "patient", "cell_type"]
)

# %%
sc.pp.normalize_total(pb_adata, target_sum=1e6)
sc.pp.log1p(pb_adata, base=2)

# %%
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    de_res_neutro = sh.compare_groups.lm.test_lm(
        pb_adata, "~ C(cell_type, Sum) + patient", groupby="cell_type", n_jobs=16
    )

# %%
cpdb_res = cpdba.significant_interactions(
    de_res_neutro.dropna(),
    max_pvalue=0.1,
    fc_col="coef",
    min_abs_fc=1,
    gene_symbol_col="variable",
    de_genes_mode="ligand",
)
cpdba.plot_result(
    cpdb_res.loc[lambda x: x["cell_type_major"].isin(["Tumor cells", "T cell CD8"])],
    fc_col="coef",
    title="Neutrophil clusters: outgoing interactions, ligand FDR < 0.1; abs(log2FC) > 1",
    aggregate=False,
    cluster="heatmap",
    de_genes_mode="ligand",
    clip_fc_at=(-3, 3),
).configure_concat(spacing=0)

# %%
cpdb_res = cpdba.significant_interactions(
    de_res_neutro.dropna(),
    max_pvalue=0.1,
    fc_col="coef",
    min_abs_fc=1,
    gene_symbol_col="variable",
    de_genes_mode="receptor",
)
cpdba.plot_result(
    cpdb_res.loc[lambda x: x["cell_type_major"].isin(["Tumor cells"])],
    fc_col="coef",
    title="Neutrophil clusters: incoming interactions, receptor FDR < 0.1; abs(log2FC) > 1",
    aggregate=False,
    cluster="heatmap",
    de_genes_mode="receptor",
    clip_fc_at=(-3, 3),
    label_limit=100,
).configure_concat(spacing=-80)

# %%

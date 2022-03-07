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
import pandas as pd
from scanpy_helpers.annotation import AnnotationHelper
import scanpy_helpers as sh
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import matplotlib
import altair as alt
from pathlib import Path
from tqdm.auto import tqdm
from nxfvars import nxfvars

# %%
ah = AnnotationHelper()

# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %%
signature_file = nxfvars.get(
    "neutro_signatures", "../../tables/gene_annotations/neutro_signatures.csv"
)
adata_file = nxfvars.get(
    "adata_in",
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
signatures = pd.read_csv(signature_file, comment="#").loc[
    lambda x: x["signature"] != "sig_pd1"
]

# %%
adata = sc.read_h5ad(adata_file)

# %% [markdown]
# # Neutrophil subset

# %%
adata_n = adata[
    (adata.obs["cell_type_coarse"] == "Neutrophils")
    & (adata.obs["condition"].isin(["LSCC", "LUAD", "NSCLC"])),
    :,
].copy()

# %%
ah.reprocess_adata_subset_scvi(
    adata_n, use_rep="X_scANVI", leiden_res=0.5, n_neighbors=15
)

# %%
sc.pl.umap(
    adata_n,
    color=[
        "leiden",
        "origin",
        "condition",
        "tumor_stage",
    ],
    wspace=0.5,
    ncols=2,
)

# %%
sc.pl.umap(
    adata_n,
    color=[
        "cell_type",
        "leiden",
        "origin",
        "condition",
        "tumor_stage",
        "sex",
        "dataset",
    ],
    wspace=0.5,
    ncols=3,
)

# %% [markdown]
# ## subclustering

# %%
sc.pl.umap(
    adata_n,
    color=["FCGR3B", "CXCL2", "OLR1", "IFIT1", "CCL3", "CDK1"],
    cmap="inferno",
    size=20,
    ncols=3,
)

# %% [markdown]
# ## Distribution of clusters across datasets and patients

# %%
patient_fracs = (
    adata_n.obs.groupby(["leiden"])
    .apply(lambda x: x["patient"].value_counts(normalize=True))
    .unstack()
    .melt(ignore_index=False, var_name="patient", value_name="fraction")
    .reset_index()
)

dataset_fracs = (
    adata_n.obs.groupby(["leiden"])
    .apply(lambda x: x["dataset"].value_counts(normalize=True))
    .unstack()
    .melt(ignore_index=False, var_name="dataset", value_name="fraction")
    .reset_index()
)

# %%
alt.Chart(dataset_fracs).mark_bar().encode(x="fraction", y="leiden", color="dataset")

# %%
alt.Chart(patient_fracs).mark_bar().encode(x="fraction", y="leiden", color="patient")

# %% [markdown]
# ## Characterization of clusters

# %%
pb_n = sh.pseudobulk.pseudobulk(adata_n, groupby=["leiden", "patient"])

# %%
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n)

# %%
pb_n = pb_n[pb_n.obs["leiden"] != "7", :].copy()

# %%
sc.tl.rank_genes_groups(pb_n, groupby="leiden", method="wilcoxon", use_raw=False)

# %%
sc.pl.rank_genes_groups_matrixplot(pb_n, dendrogram=False)

# %%
sc.pl.umap(adata_n, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
marker_genes = pd.DataFrame(pb_n.uns["rank_genes_groups"]["names"])

# %%
# markers from zillionis/klein
sc.pl.umap(
    adata_n, color=["IFIT1", "MMP8", "PALD1", "CXCL3", "CTSC", "CCL3"], cmap="inferno"
)

# %%
for i, cl in enumerate(pb_n.obs["leiden"].cat.categories):
    print(f"Leiden: {cl}")
    sc.pl.umap(
        adata_n,
        color=marker_genes.iloc[:5, i].tolist(),
        cmap="inferno",
        size=20,
        ncols=5,
    )

# %% [markdown]
# ## Annotate clusters

# %%
ah.annotate_cell_types(
    adata_n,
    {
        "TAN CTSC": [1],  # antigen presentation; also in zillionis/klein
        "NAN S100A12": [4],  # also in zillionis/klein
        "NAN CXCR2": [0],
        "TAN IFIT": [6],  # also in zillionis/klein
        "TAN RPL5": [3],
        "TAN CXCL2": [5, 2],
    },
)

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(
        adata_n,
        color="cell_type",
        legend_loc="on data",
        legend_fontoutline=2,
        frameon=False,
        legend_fontsize=10,
        size=40,
    )

# %% [markdown]
# ## Heatmaps by neutrophil clusters

# %%
sc.tl.rank_genes_groups(adata_n, "cell_type")

# %%
sc.pl.umap(
    adata_n,
    color=["CXCR2", "S100A12", "CTSC", "CXCL2", "IFIT1", "RPL5"],
    cmap="inferno",
    size=20,
    ncols=3,
    frameon=False,
)

# %%
sc.pl.dotplot(
    adata_n,
    var_names=["CXCR2", "S100A12", "CTSC", "CXCL2", "IFIT1", "RPL5"],
    groupby="cell_type",
)

# %%
sc.pl.rank_genes_groups_dotplot(adata_n, dendrogram=False)

# %%
pb_n = sh.pseudobulk.pseudobulk(adata_n, groupby=["cell_type", "patient"])

# %%
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n)

# %%
pb_n = pb_n[pb_n.obs["cell_type"] != "other", :].copy()

# %%
sc.tl.rank_genes_groups(pb_n, groupby="cell_type", method="wilcoxon", use_raw=False)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=["CXCR2", "S100A12", "CTSC", "CXCL2", "IFIT1", "RPL5"],
    groupby="cell_type",
    cmap="bwr",
)

# %%
sc.pl.rank_genes_groups_matrixplot(pb_n, dendrogram=False, cmap="bwr")

# %%
signature_dict = {}
for sig in signatures["signature"].unique():
    signature_dict[sig] = signatures.loc[
        lambda x: x["signature"] == sig, "gene_symbol"
    ].tolist()

# %%
sc.pl.matrixplot(pb_n, var_names=signature_dict, cmap="bwr", groupby="cell_type")

# %%
for sig, genes in signature_dict.items():
    sc.tl.score_genes(pb_n, genes, score_name=sig)
    sc.tl.score_genes(adata_n, genes, score_name=sig)

# %%
sc.pl.matrixplot(
    pb_n, var_names=list(signature_dict.keys()), cmap="bwr", groupby="cell_type"
)

# %%
sc.pl.matrixplot(
    adata_n[adata_n.obs["cell_type"] != "other"],
    var_names=list(signature_dict.keys()),
    cmap="bwr",
    groupby="cell_type",
)

# %%
sc.pl.umap(adata_n, color=list(signature_dict.keys()), cmap="inferno", size=20, ncols=3)

# %% [markdown]
# ## Merge annotation into atlas

# %%
adata_n.obs["cell_type_neutro"] = adata_n.obs["cell_type"]
adata.obs["cell_type_neutro"] = adata.obs["cell_type_major"]
ah.integrate_back(adata, adata_n, variable="cell_type_neutro")

# %% [markdown]
# ## UKIM-V datasets only

# %%
adata_n_ukimv = adata_n[adata_n.obs["dataset"].str.startswith("UKIM-V"), :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_n_ukimv, use_rep="X_scANVI")

# %%
adata_n_ukimv.obs

# %%
with plt.rc_context({"figure.figsize": (4, 4), "figure.dpi": 300}):
    sc.pl.umap(
        adata_n,
        color="origin",
        frameon=False,
        size=20,
        groups=["normal_adjacent", "tumor_metastasis", "tumor_primary"],
    )
    sc.pl.umap(adata_n, color="condition", frameon=False, size=20)

# %%
tmp_df = (
    adata_n.obs.groupby(["cell_type"])
    .apply(lambda x: x["origin"].value_counts(normalize=True))
    .unstack()
    .melt(var_name="origin", value_name="fraction", ignore_index=False)
    .reset_index()
    .loc[lambda x: x["cell_type"] != "other"]
)

# %%
alt.Chart(tmp_df).mark_bar().encode(x="fraction", y="cell_type", color="origin")

# %%
tmp_df = (
    adata_n.obs.groupby(["cell_type"])
    .apply(lambda x: x["condition"].value_counts(normalize=True))
    .unstack()
    .melt(var_name="condition", value_name="fraction", ignore_index=False)
    .reset_index()
    .loc[lambda x: x["cell_type"] != "other"]
)

# %%
alt.Chart(tmp_df).mark_bar().encode(x="fraction", y="cell_type", color="condition")

# %% [markdown]
# # Save result

# %%
adata_n.write_h5ad(f"{artifact_dir}/adata_neutrophil_clusters.h5ad")
adata.write_h5ad(f"{artifact_dir}/full_atlas_neutrophil_clusters.h5ad")

# %%

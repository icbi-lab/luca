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
sc.pl.umap(adata_n, color="leiden", legend_loc="on data", legend_fontoutline=2)

# %%
ah.annotate_cell_types(
    adata_n,
    {
        "TAN-2": [4],  # antigen presentation; also in zillionis/klein
        "NAN-1": [2],  # also in zillionis/klein
        "NAN-2": [1],
        "TAN-4": [5],  # also in zillionis/klein
        "TAN-1": [3, 6],
        "TAN-3": [0],
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
# ## Merge annotation into atlas

# %%
adata_n.obs["cell_type_neutro"] = adata_n.obs["cell_type"]
adata.obs["cell_type_neutro"] = adata.obs["cell_type_major"]
ah.integrate_back(adata, adata_n, variable="cell_type_neutro")

# %%
adata.obs["cell_type_neutro_coarse"] = adata.obs["cell_type_major"].astype(str)
adata.obs.loc[
    adata.obs["cell_type_neutro"].str.startswith("NAN"), "cell_type_neutro_coarse"
] = "NAN"
adata.obs.loc[
    adata.obs["cell_type_neutro"].str.startswith("TAN"), "cell_type_neutro_coarse"
] = "TAN"

# %% [markdown]
# # Save result

# %%
adata_n.write_h5ad(f"{artifact_dir}/adata_neutrophil_clusters.h5ad")
adata.write_h5ad(f"{artifact_dir}/full_atlas_neutrophil_clusters.h5ad")

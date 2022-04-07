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
adata_file = nxfvars.get(
    "adata_in",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
neutro_clustering = nxfvars.get("neutro_clustering", "../../tables/neutrophil_clustering.csv")
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
adata = sc.read_h5ad(adata_file)

# %%
neutrophil_clustering = pd.read_csv(neutro_clustering, index_col=0)

# %% [markdown]
# # Neutrophil subset

# %%
adata_n = adata[
    (adata.obs["cell_type_coarse"] == "Neutrophils")
    & (adata.obs["condition"].isin(["LUSC", "LUAD", "NSCLC NOS"])),
    :,
].copy()

# %%
adata_n.shape

# %%
ah.reprocess_adata_subset_scvi(
    adata_n, use_rep="X_scANVI", leiden_res=0.5, n_neighbors=20
)

# %%
# We use a predefined clustering for neutrophils here. 
# 
# It was generated using exactly this notebook, but using an earlier version of the atlas. 
# After fixing an unrelated mistake in the atlas metadata, the scANVI embedding has changed slightly. 
# While the clustering is qualitativly similar, the clusters don't match 1:1 with the previous version
# and at this point, we didn't want to change downstream results anymore. 
adata_n.obs["leiden"] = neutrophil_clustering["leiden"].astype(str)
adata_n = adata_n[~adata_n.obs["leiden"].isnull(), :].copy()

# %%
adata_n.shape

# %%
# flip UMAP y axis to be visually consistent with previous iterations of the dataset
adata_n.obsm["X_umap"][:, 1] = np.max(adata_n.obsm["X_umap"][:, 1]) - adata_n.obsm["X_umap"][:, 1]

# %%
sc.pl.umap(adata_n, color="dataset")

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
        "NAN-1": [4],  # also in zillionis/klein, S100A12
        "NAN-2": [0],
        "TAN-1": [3], 
        "TAN-2": [2], 
        "TAN-3": [1],  
        "TAN-4": [5], 
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
# Get rid of the ~10 cells that don't have a neutrophil cluster assigned
adata = adata[adata.obs["cell_type_neutro"] != "Neutrophils", :]

# %%
sc.pl.umap(adata, color="cell_type_neutro")

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

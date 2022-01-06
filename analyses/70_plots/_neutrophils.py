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
import scanpy as sc
import pandas as pd
from scanpy_helpers.annotation import AnnotationHelper

# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
adata = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad"
)

# %%
adata_t = adata[adata.obs["cell_type_coarse"] == "Tumor cells", :].copy()
adata_n = adata[adata.obs["cell_type_coarse"] == "Neutrophils", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_n, use_rep="X_scANVI")

# %%
sc.pl.umap(adata_n, color=["cell_type",  "origin", "condition", "tumor_stage", "sex", "dataset"], wspace=0.5)

# %% [markdown]
# ### Compare cytosig and progeny
#
# ### Compare DE genes using pair plots

# %% [markdown]
# ---

# %%
adata.obs.groupby(["dataset", "patient", "sample", "origin"], observed=True).size().reset_index(
    name="n_cells"
).sort_values("n_cells", ascending=False).to_csv("./granulocyte_count_per_patient.tsv", sep="\t")

# %%
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)

# %%
sc.tl.pca(adata, use_highly_variable=True)

# %%
sc.pp.neighbors(adata)

# %%
sc.tl.umap(adata)

# %%
sc.pl.umap(adata, color=["cell_type", "dataset"])

# %%
adata

# %%

# %% [markdown] tags=[]
# ---
#
# # Neutrophils

# %%
adata_neutro = adata[adata.obs["cell_type"] == "Granulocytes", :].copy()

# %%
sc.pp.neighbors(adata_neutro, use_rep="X_scANVI")
sc.tl.leiden(adata_neutro, resolution=0.5)
sc.tl.umap(adata_neutro)

# %%
sc.pl.umap(adata_neutro, color=["origin", "leiden", "dataset"], wspace=0.5)

# %%
sc.pl.umap(
    adata_neutro,
    color=adata.obs.columns[adata.obs.columns.str.startswith("scissor")],
    wspace=0.5,
    ncols=3,
)

# %%
neutro_ukimv = adata_neutro[adata_neutro.obs["dataset"] == "UKIM-V", :]

# %%
sc.pl.umap(
    neutro_ukimv,
    color=["origin", "leiden", "patient", "dataset", "VEGFA"],
    wspace=0.5,
    ncols=2,
)

# %%
sc.pl.dotplot(
    neutro_ukimv,
    groupby=["patient", "origin"],
    var_names="VEGFA",
    title="tumor_primary",
    vmin=0, 
    vmax=1
)

# %%
sc.pl.dotplot(
    neutro_ukimv[neutro_ukimv.obs["origin"] == "normal_adjacent", :],
    groupby="patient",
    var_names="VEGFA",
    title="normal_adjacent",
    vmin=0, 
    vmax=1
)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata_neutro, color=["VEGFA"], cmap="inferno")

# %%
sc.pl.dotplot(
    adata,
    groupby=["patient", "origin"],
    var_names="VEGFA",
    title="tumor_primary",
    vmin=0, 
    vmax=1
)

# %%
adata_neutro.obs.groupby(
    ["dataset", "patient", "origin"], observed=True
).size().reset_index(name="n")

# %%
sc.pl.umap(adata_neutro, color=["origin", "dataset"])

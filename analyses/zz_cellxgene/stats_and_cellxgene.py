# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
import scanpy as sc
from nxfvars import nxfvars

# %%
adata_in = nxfvars.get(
    "adata_in",
    "../../data/20_integrate_scrnaseq_data/29_annotate_cell_types_fine/artifacts/adata_annotated_fine.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/data/scratch/sturm/tmp/")

# %% [markdown]
# ### Prepare datasets for cellxgene
#
# -> all genes in `X`, normalized values in `X`. 

# %%
adata = sc.read_h5ad(adata_in)

# %% [markdown]
# ### Dimensions (number of cells, number of genes)

# %%
adata.raw.shape

# %% [markdown]
# ### Overview of celltypes

# %%
sc.set_figure_params(figsize=(10, 10))
sc.pl.umap(
    adata,
    color=["cell_type_coarse", "cell_type"],
    legend_loc="on data",
    legend_fontoutline=2,
    legend_fontsize=8,
    size=1,
)

# %%
sc.set_figure_params(figsize=(6, 6))
sc.pl.umap(adata, color=["dataset", "condition", "origin"], wspace=0.5)

# %% [markdown]
# ### count by cell-type

# %%
adata.obs["cell_type_coarse"].value_counts()

# %%
adata.obs["cell_type"].value_counts()

# %% [markdown]
# ### count by origin

# %%
adata.obs["origin"].value_counts()

# %% [markdown]
# ### count by dataset

# %%
adata.obs["dataset"].value_counts()

# %%
print(f"total datasets: {adata.obs['dataset'].nunique()}")

# %% [markdown]
# ### count by condition

# %%
adata.obs["condition"].value_counts()

# %% [markdown]
# ### count by tissue

# %%
adata.obs["tissue"].value_counts()

# %% [markdown]
# ### number of patients

# %%
print(f"total patients: {adata.obs['patient'].nunique()}")

# %% [markdown]
# ### number of patients per origin

# %%
adata.obs.loc[:, ["origin", "patient"]].drop_duplicates()["origin"].value_counts()

# %% [markdown]
# ### Export for cellxgene

# %%
adata_cellxgene = sc.AnnData(
    var=adata.raw.var, X=adata.raw.X, obs=adata.obs, obsm=adata.obsm
)

# %%
adata_cellxgene.write_h5ad(f"{artifact_dir}/adata_cellxgene.h5ad")

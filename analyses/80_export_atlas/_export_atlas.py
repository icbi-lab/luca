# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:conda-2020-pircher-cellxgene-schema]
#     language: python
#     name: conda-env-conda-2020-pircher-cellxgene-schema-py
# ---

# %%
import scanpy as sc
from nxfvars import nxfvars
import pandas as pd

# %% [markdown]
# Export final atlas for sharing. 
#  * Use ontologies
#  * Remove unnecessary data
#  * make sure it passes the cellxgene schema
#  * output documentation

# %%
path_core_atlas = nxfvars.get(
    "core_atlas",
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)
path_extended_atlas = nxfvars.get(
    "extended_atlas",
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad"
)
path_neutrophil_atlas = nxfvars.get(
    "neutrophil_atlas",
    "../../data/30_downstream_analyses/neutrophils/subclustering/artifacts/full_atlas_neutrophil_clusters.h5ad")

# %%
core_atlas = sc.read_h5ad(path_core_atlas)

# %%
extended_atlas = sc.read_h5ad(path_extended_atlas)

# %%
neutrophil_atlas = sc.read_h5ad(path_neutrophil_atlas)

# %%
core_atlas.shape

# %%
extended_atlas.shape

# %%
neutrophil_atlas.shape

# %%
pd.set_option("display.max_columns", None)

# %% [markdown]
# ## Merge all information, clean up object

# %% [markdown]
# ### Add Neutrophil annotations

# %%
neutrophil_atlas.obs["cell_type_neutro_coarse"]

# %%
# TODO currently not the case, but should be in the final version
assert extended_atlas.shape == neutrophil_atlas.shape

# TODO set with copy warning
for c in ["cell_type_neutro", "cell_type_neutro_coarse"]:
    extended_atlas.obs[c] = neutrophil_atlas.obs[c]

# %% [markdown]
# ### Remove unnecessary columns

# %%
extended_atlas.obs = extended_atlas.obs.loc[:, lambda x: ~x.columns.str.startswith("_")]

# %% [markdown]
# ## Reannotate with ontologies

# %% [markdown]
# ## Export h5ads

# %% [markdown]
# ## Validate h5ads

# %%
# TODO
# !cellxgene-schema validate {path_core_atlas}

# %%

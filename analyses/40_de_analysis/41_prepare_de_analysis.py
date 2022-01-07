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

# %% [markdown]
# ## Split anndata into subsets that only contain the required data for the DE analysis. 

# %%
import scanpy as sc
import pandas as pd
from nxfvars import nxfvars

# %%
# only contains 6000 most DE genes
adata = sc.read_h5ad(
    nxfvars.get(
        "input_adata",
        "../../data/20_build_atlas//annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
    )
)

artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
sc.set_figure_params(figsize=(9, 9))

# %%
sc.pl.umap(adata, color=["cell_type"], size=1)

# %%
obs = adata.obs.loc[
    :, ["sample", "patient", "tissue", "origin", "condition", "dataset", "cell_type"]
]

# %%
obs

# %%
pd.set_option("display.max_rows", 1000)

# %% [markdown]
# # Comparison 1: paired tumor vs normal

# %%
patients_w_normal = set(
    obs.loc[obs["origin"].isin(["normal_adjacent", "normal"]), :]["patient"]
)

# %%
patients_w_primary_tumor = set(obs.loc[obs["origin"] == "tumor_primary", :]["patient"])

# %%
len(patients_w_normal)

# %%
len(patients_w_primary_tumor)

# %%
len(patients_w_normal & patients_w_primary_tumor)

# %%
obs.groupby(
    ["patient", "tissue", "origin", "condition"], observed=True
).size().reset_index(name="count").sort_values(["patient", "origin"]).loc[
    lambda x: x["patient"].isin(patients_w_normal & patients_w_primary_tumor), :
]

# %%
adata_tumor_normal = adata[
    adata.obs["patient"].isin(patients_w_normal & patients_w_primary_tumor)
    & adata.obs["origin"].isin(["normal_adjacent", "normal", "tumor_primary"]),
    :,
]

# %%
adata_tumor_normal

# %%
sc.pl.umap(adata_tumor_normal, color="cell_type")

# %% [markdown]
# # Comparison 2: LUAD vs LSCC
#
# primary tumor samples from LUAD and LSCC patients

# %%
adata_luad_lusc = adata[
    adata.obs["condition"].isin(["LUAD", "LSCC"])
    & (adata.obs["origin"] == "tumor_primary"),
    :,
].copy()

# %%
sc.pl.umap(adata_luad_lusc, color="cell_type")

# %% [markdown]
# # Write output

# %%
adata_tumor_normal.write_h5ad(f"{artifact_dir}/adata_tumor_normal.h5ad")

# %%
adata_luad_lusc.write_h5ad(f"{artifact_dir}/adata_luad_lscc.h5ad")

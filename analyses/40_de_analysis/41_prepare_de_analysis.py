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
        "../../data/30_downstream_analyses/03_update_annotation/artifacts/full_atlas_merged.h5ad",
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
# ## Ensure patients are only in one dataset
#
# For the leader/merad datasets, some patients have beend sequenced across multiple batches.
# We exclude all duplicate patients from other batches than the main batch ("10x_3p_v2").
# The same patient across multiple datasets would require more advanced models, which are not worth it for just 4 patients.

# %%
patients_in_multiple_datasets = (
    adata.obs.groupby("patient")
    .apply(lambda x: x["dataset"].nunique())
    .sort_values(ascending=False)[lambda x: x > 1]
)

# %%
patients_in_multiple_datasets

# %%
adata.obs.loc[
    lambda x: x["patient"].isin(patients_in_multiple_datasets.index),
    ["patient", "sample", "dataset"],
].drop_duplicates()

# %%
adata = adata[
    (adata.obs["dataset"] != "Leader_Merad_2021_10x_3p_v2_digest-deadcell_cite")
    & ~(
        (adata.obs["dataset"] == "Leader_Merad_2021_10x_3p_v2_beads_cite")
        & (adata.obs["patient"] == "Leader_Merad_2021_581")
    )
    & ~(
        (adata.obs["dataset"] == "Leader_Merad_2021_10x_5p_v1_beads")
        & (adata.obs["patient"] == "Leader_Merad_2021_522")
    ),
    :,
].copy()

# %%
assert all(adata.obs.groupby("patient").apply(lambda x: x["dataset"].nunique()) == 1)

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
# # Comparison 2: LUAD vs LUSC
#
# primary tumor samples from LUAD and LUSC patients

# %%
adata_luad_lusc = adata[
    (adata.obs["origin"] == "tumor_priLUSC"),
    :,
].copy()

# %%
sc.pl.umap(adata_luad_lusc, color="cell_type")

# %% [markdown]
# # Write output

# %%
adata_tumor_normal.write_h5ad(f"{artifact_dir}/adata_tumor_normal.h5ad")

# %%
adata_luad_lusc.write_h5ad(f"{artifact_dir}/adata_primary_tumor.h5ad")

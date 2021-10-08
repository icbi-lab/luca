# ---
# jupyter:
#   jupytext:
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
        "../../data/20_integrate_scrnaseq_data/29_annotate_cell_types_fine/artifacts/adata_annotated_fine.h5ad",
    )
)

artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")

# %%
sc.set_figure_params(figsize=(9, 9))

# %%
sc.pl.umap(adata, color=["cell_type"], size=3)

# %%
obs = adata.obs.loc[
    :, ["sample", "patient", "tissue", "origin", "condition", "dataset", "cell_type"]
]

# %%
obs

# %%
pd.set_option("display.max_rows", 1000)

# %% [markdown]
# ### Comparison 1: paired tumor vs normal

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

# %%
adata_tumor_normal.write_h5ad(f"{artifact_dir}/adata_tumor_normal.h5ad")

# %% [markdown]
# ### Comparison 2: Healthy control vs. Tumor; COPD vs Tumor; COPD vs. healthy control

# %%
adata.obs["condition"].unique()

# %%
obs.loc[
    obs["condition"].isin(["NSCLC", "LSCC", "LUAD"]), "dataset"
].drop_duplicates().tolist()

# %%
obs.loc[obs["condition"] == "COPD", "dataset"].drop_duplicates().tolist()

# %%
obs.loc[obs["condition"] == "healthy_control", "dataset"].drop_duplicates().tolist()

# %%
adata.obs["condition"] = [
    "NSCLC" if c in ["LSCC", "LUAD"] else c for c in adata.obs["condition"]
]

# %%
adata.obs["condition_origin"] = [
    f"{condition}_{origin}"
    for condition, origin in zip(adata.obs["condition"], adata.obs["origin"])
]

# %%
adata_healthy_copd_tumor = adata[
    adata.obs["condition_origin"].isin(
        ["COPD_normal", "NSCLC_tumor_primary", "healthy_control_normal"]
    ), :
]

# %%
adata_healthy_copd_tumor

# %%
adata_healthy_copd_tumor.obs["condition_origin"].unique()

# %%
sc.pl.umap(adata_healthy_copd_tumor, color="cell_type")

# %%
adata_healthy_copd_tumor.write_h5ad(f"{artifact_dir}/adata_healthy_copd_nsclc.h5ad")

# %%

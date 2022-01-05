# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
from nxfvars import nxfvars
from scanpy_helpers.annotation import AnnotationHelper
from datetime import datetime
import pandas as pd
import numpy as np

# %%
sc.settings.set_figure_params(figsize=(6, 6))

# %%
path_adata_fine = nxfvars.get(
    "adata_annotated_fine",
    "../../data/20_build_atlas/annotate_datasets/32_cell_types_fine/artifacts/adata_annotated_fine.h5ad",
)
path_adata_epi = nxfvars.get(
    "adata_epi",
    "../../data/20_build_atlas/annotate_datasets/33_cell_types_epi/artifacts/adata_epithelial.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/data/scratch/sturm/tmp/")

# %%
ah = AnnotationHelper()

# %%
adata_epi = sc.read_h5ad(path_adata_epi)
adata = sc.read_h5ad(path_adata_fine)

# %% [markdown]
# # Merge cell-type annotations

# %%
sc.pl.umap(adata, color="cell_type")

# %%
sc.pl.umap(adata_epi, color="cell_type")

# %%
ah.integrate_back(adata, adata_epi)

# %%
adata.obs["cell_type_tumor"] = adata.obs["cell_type"]
adata.obs["cell_type"] = [
    "Tumor cells" if x.startswith("Tumor") else x for x in adata.obs["cell_type_tumor"]
]

# %%
adata_epi.obs["cell_type_tumor"] = adata.obs["cell_type_tumor"]
adata_epi.obs["cell_type"] = adata.obs["cell_type"]

# %%
adata = adata[
    ~adata.obs["cell_type"].isin(["Neuronal cells", "Hepatocytes", "Hemoglobin+"]), :
]

# %%
adata.obs.columns

# %%
for col in adata.obs.columns:
    if col.startswith("_"):
        del adata.obs[col]

del adata.obs["leiden_1.00"]

# %%
adata.obs.columns

# %% [markdown]
# # Add missing dataset metadata

# %%
wu = pd.read_excel(
    "../../tables/additional_patient_metadata/wu_et_al_additional_patient_metadata.xlsx"
)

# %%
wu2 = (
    wu.loc[:, ["Patient", "FINAL"]]
    .rename(columns={"Patient": "patient", "FINAL": "pathology"})
    .assign(patient=lambda x: [f"Wu_Zhou_2021_NSCLC_{p}" for p in x["patient"]])
)

# %%
mayn = pd.read_csv(
    "../../tables/additional_patient_metadata/cell_metadata_maynard.csv.gz"
)

# %%
pd.set_option("display.max_columns", 200)

# %%
mayn2 = (
    mayn.loc[
        :,
        ["patient_id", "gender", "smokingHx", "histolgy", "stage.at.dx", "driver_gene"],
    ]
    .drop_duplicates(ignore_index=True)
    .rename(
        columns={
            "patient_id": "patient",
            "gender": "sex",
            "smokingHx": "ever_smoker",
            "histolgy": "pathology",
            "stage.at.dx": "uicc_stage",
            "driver_gene": "driver_genes",
        }
    )
    .assign(
        ever_smoker=lambda x: [
            {"Never": "no", "Former": "yes", "Current": "yes"}[s]
            for s in x["ever_smoker"]
        ],
        patient=lambda x: [f"Maynard_Bivona_2020_NSCLC_{p}" for p in x["patient"]],
    )
)

# %%
add_meta = (
    pd.concat(
        [
            pd.read_excel(
                "../../tables/additional_patient_metadata/additional_metadata.xlsx"
            ),
            mayn2,
            wu2,
        ]
    )
    .drop_duplicates()
    .rename(columns={"pathology": "condition"})
    .assign(
        condition= lambda x: [
            {
                "control": "healthy_control",
                "copd": "COPD",
                "squamous": "LSCC",
                "large cell": "LCLC",
                "adeno": "LUAD",
                "pleomorphic": "PPC",
                "adenocarcinoma": "LUAD",
                "lusc": "LSCC",
                "luad": "LUAD",
                "nsclc": "NSCLC",
            }[c.lower()] if not pd.isnull(c) else None for c in x["condition"]
        ],
        sex = lambda x: x.sex.str.lower().str.strip()
    )
)

# %%
add_meta.sex.unique()

# %%
add_meta.condition.unique()

# %%
add_meta

# %%
add_meta.loc[~add_meta["patient"].isin(adata.obs["patient"]), :]

# %%
meta_from_adata = adata.obs.loc[:, ["patient", "sex", "condition"]].drop_duplicates(
    ignore_index=True
)

# %%
meta_from_adata

# %%
pd.set_option("display.max_rows", 200)

# %%
meta_from_adata.set_index("patient").join(
    add_meta.set_index("patient"), how="outer", lsuffix="_atlas", rsuffix="_new"
).sort_index(axis="columns").sort_index().to_excel("../../tables/additional_patient_metadata/patient_metadata.xlsx")

# %%
meta_from_adata.loc[~meta_from_adata["patient"].isin(add_meta["patient"]), :]

# %% [markdown]
# # Dataset overview

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
sc.pl.umap(adata, color=["condition", "origin", "dataset"], wspace=0.3)

# %% [markdown]
# ### Dimensions (number of cells, number of genes)

# %%
adata.raw.shape

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
# ## Export data

# %% [markdown]
# ### Prepare datasets for cellxgene
#
# -> all genes in `X`, normalized values in `X`. 

# %%
date_now = datetime.now().date().isoformat()

# %%
adata_cellxgene = sc.AnnData(
    var=adata.raw.var, X=adata.raw.X, obs=adata.obs, obsm=adata.obsm
)

# %%
adata_cellxgene_epi = sc.AnnData(
    var=adata_epi.raw.var, X=adata_epi.raw.X, obs=adata_epi.obs, obsm=adata_epi.obsm
)

# %%
adata.write_h5ad(f"{artifact_dir}/full_atlas_annotated.h5ad")
adata_epi.write_h5ad(f"{artifact_dir}/epithelial_cells_annotated.h5ad")

# %%
adata_cellxgene.write_h5ad(f"{artifact_dir}/full_atlas_{date_now}.h5ad")
adata_cellxgene_epi.write_h5ad(f"{artifact_dir}/epithelial_cells_{date_now}.h5ad")

# %%

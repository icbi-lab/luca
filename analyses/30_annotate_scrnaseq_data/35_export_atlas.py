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
path_patient_metadata = nxfvars.get(
    "patient_metadata",
    "../../tables/additional_patient_metadata/patient_metadata_corrected.xlsx"
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
].copy()

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
additional_patient_metadata = pd.read_excel(
    path_patient_metadata
).assign(ever_smoker=lambda x: x["ever_smoker"].str.strip().str.lower())
additional_patient_metadata

# %%
for col in additional_patient_metadata.columns[2:]:
    print(additional_patient_metadata[col].unique().tolist())

# %%
additional_metadata_by_cell = (
    adata.obs.drop(["condition", "sex"], axis="columns")
    .reset_index()
    .merge(additional_patient_metadata, on="patient")
    .set_index("index")
)

# %%
adata.obs.drop(["condition", "sex"], axis="columns", inplace=True)
for col in additional_patient_metadata.columns:
    if col != "patient":
        adata.obs.insert(1, col, additional_metadata_by_cell[col])

# %%
adata.obs["tumor_stage"] = np.nan
adata.obs.loc[adata.obs["uicc_stage"].isin(["I", "II", "IA"]), "tumor_stage"] = "early"
adata.obs.loc[adata.obs["uicc_stage"].isin(["III", "IV", "IIIA", "IIIB", "III or IV"]), "tumor_stage"] = "late"
assert np.all(pd.isnull(adata.obs["uicc_stage"]) == pd.isnull(adata.obs["tumor_stage"]))

# %% [markdown]
# ### Simplify dataset names

# %%
remap_datasets = {d: "_".join(d.split("_")[:-1]) for d in adata.obs["dataset"].unique()}
remap_datasets["UKIM-V"] = "UKIM-V"
remap_datasets["Travaglini_Krasnow_2020_Lung_SS2"] = "Travaglini_Krasnow_2020"
remap_datasets["Lambrechts_2018_LUAD_6149v1"] = "Lambrechts_Thienpont_2018_6149v1"
remap_datasets["Lambrechts_2018_LUAD_6149v2"] = "Lambrechts_Thienpont_2018_6149v2"
remap_datasets["Lambrechts_2018_LUAD_6653"] = "Lambrechts_Thienpont_2018_6653"

# %%
pd.DataFrame().from_dict(remap_datasets, orient="index", columns=["new"])

# %%
adata.obs.drop("batch", axis="columns", inplace=True)

# %%
pd.set_option("display.max_columns", 200)

# %%
for col in ["sample", "patient", "dataset"]:
    for old, new in remap_datasets.items():
        adata.obs[col] = adata.obs[col].str.replace(old, new, regex=False)

# %%
adata.obs

# %% [markdown]
# # Adjust cell-type annotations

# %%
adata.obs["cell_type_coarse"] = (
    adata.obs["cell_type_coarse"].str.replace("Granulocytes", "Neutrophils")
).astype(str)
adata.obs["cell_type"] = adata.obs["cell_type"].str.replace(
    "Granulocytes", "Neutrophils"
)

# %%
adata.obs.loc[adata.obs["cell_type"] == "NK cell", "cell_type_coarse"] = "NK cell"

# %%
adata.obs.loc[
    adata.obs["cell_type"].isin(
        ["Macrophage FABP4+", "Macrophage", "Monocyte", "myeloid dividing"]
    ),
    "cell_type_coarse",
] = "Macrophage/Monocyte"
adata.obs.loc[
    adata.obs["cell_type"].isin(["cDC1", "cDC2", "DC mature"]), "cell_type_coarse"
] = "cDC"

# %% [markdown]
# ### cell_type_major
#
# an annotation that is more fine-grained than coarse, more coarse grained than fine, and ignore some special categories, such as dividing cells. 
# This may be useful for some downstream analyses. 

# %%
adata.obs["cell_type_major"] = adata.obs["cell_type"].astype(str)
adata.obs.loc[
    adata.obs["cell_type_coarse"] == "Stromal",
    "cell_type_major",
] = "Stromal"
adata.obs.loc[
    adata.obs["cell_type_coarse"] == "Endothelial cell",
    "cell_type_major",
] = "Endothelial cell"
adata.obs.loc[
    adata.obs["cell_type_major"].isin(
        [
            "other (T assoc.)",
            "T cell dividing",
            "B cell dividing",
            "myeloid dividing",
            "ROS1+ healthy epithelial",
        ]
    ),
    "cell_type_major",
] = "other"

# %% [markdown]
# # Dataset overview

# %% [markdown]
# ### Overview of celltypes

# %%
sc.set_figure_params(figsize=(10, 10))
sc.pl.umap(
    adata,
    color=["cell_type_coarse", "cell_type", "cell_type_major"],
    legend_loc="on data",
    legend_fontoutline=2,
    legend_fontsize=8,
    size=1,
)

# %%
sc.set_figure_params(figsize=(6, 6))
sc.pl.umap(
    adata,
    color=["condition", "origin", "dataset", "sex", "ever_smoker", "uicc_stage", "tumor_stage"],
    wspace=0.3,
    ncols=3,
)

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

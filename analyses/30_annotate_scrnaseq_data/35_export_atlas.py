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
import warnings
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
    "../../tables/additional_patient_metadata/patient_metadata_corrected.xlsx",
)
path_platform_metadata = nxfvars.get(
    "platform_metadata",
    "../../tables/additional_patient_metadata/sequencing_platforms.csv",
)
artifact_dir = nxfvars.get(
    "artifact_dir",
    "/home/sturm/Downloads",
)

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
# exclude doublets and cell-type from distant metastases
adata = adata[
    ~(
        adata.obs["cell_type"].isin(
            ["Neuronal cells", "Hepatocytes", "Hemoglobin+", "unknown (plasma assoc.)"]
        )
        | adata.obs["cell_type"].str.contains("empty droplet")
        | adata.obs["cell_type"].str.contains("doublet")
    ),
    :,
].copy()

# %%
adata.obs["cell_type"].unique().tolist()

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
additional_patient_metadata = pd.read_excel(path_patient_metadata).assign(
    ever_smoker=lambda x: x["ever_smoker"].str.strip().str.lower()
)
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
adata.obs.loc[adata.obs["uicc_stage"].isin(["I", "II"]), "tumor_stage"] = "early"
adata.obs.loc[
    adata.obs["uicc_stage"].isin(["III", "IV", "III or IV"]),
    "tumor_stage",
] = "advanced"
assert np.all(pd.isnull(adata.obs["uicc_stage"]) == pd.isnull(adata.obs["tumor_stage"]))

# %% [markdown]
# ## Correct metadata
#
#  * fix metadata that is incorrectly assigned in the original datasets or
#  * make metadata more accessible

# %%
# rename "healthy_control" to "non-cancer"
adata.obs.loc[lambda x: x["condition"] == "healthy_control", "condition"] = "non-cancer"
# rename LSCC -> LUSC
adata.obs.loc[lambda x: x["condition"] == "LSCC", "condition"] = "LUSC"
# rename LCLC, PPC -> NOS
adata.obs.loc[
    lambda x: x["condition"].isin(["LCLC", "PPC", "NSCLC"]), "condition"
] = "NSCLC NOS"
adata.obs["condition"].value_counts()

# %%
# make dedicated column for each relevant driver mutations.
# ignoring driver mutations that are only specified in few patients
for gene in ["EGFR", "TP53", "ALK", "BRAF", "ERBB2", "KRAS", "ROS"]:
    adata.obs[f"{gene}_mutation"] = [
        (("mutated" if gene in x else "not mutated") if not pd.isnull(x) else np.nan)
        for x in adata.obs["driver_genes"]
    ]

# %%
# He_Fan has tumor and normal samples switched
adata.obs["origin"] = [
    {"tumor_primary": "normal_adjacent", "normal_adjacent": "tumor_primary"}[x]
    if dataset == "He_Fan_2021_LUAD"
    else x
    for x, dataset in zip(adata.obs["origin"], adata.obs["dataset"])
]
adata.obs["sample"] = [
    x.replace("LUAD_N", "LUAD_tmp").replace("LUAD_LUAD", "LUAD_N").replace("LUAD_tmp", "LUAD_LUAD")
    if dataset == "He_Fan_2021_LUAD"
    else x
    for x, dataset in zip(adata.obs["sample"], adata.obs["dataset"])
]

# aggregate additonal tumor origins from Lambrechts under "tumor_primary". Keep original annotation in separate column
adata.obs["origin_fine"] = adata.obs["origin"].copy()
adata.obs.loc[lambda x: x["origin"] == "tumor_middle", "origin"] = "tumor_primary"
adata.obs.loc[lambda x: x["origin"] == "tumor_edge", "origin"] = "tumor_primary"

# aggregate pleura tissue and effusion origin into pleura/effusion tissue
adata.obs["tissue"] = adata.obs["tissue"].astype(str)
adata.obs.loc[
    lambda x: (x["origin"] == "effusion") | (x["tissue"] == "pleura"), "tissue"
] = "pleura/effusion"

# categorize effusion samples as tumor metastases as it all comes from pleura metastases
adata.obs.loc[lambda x: x["origin"] == "effusion", "origin"] = "tumor_metastasis"

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    adata._sanitize()

# %%
print(adata.obs["origin"].value_counts())
print()
print(adata.obs["tissue"].value_counts())

# %%
# Normal samples from Maier_Merad should be normal_adjacent
adata.obs.loc[
    lambda x: x["dataset"].str.contains("Maier_Merad") & (x["origin"] == "normal"),
    "origin",
] = "normal_adjacent"
# Normal samples from Mayr_Schiller should be normal adjacent
adata.obs.loc[
    lambda x: x["dataset"].str.contains("Mayr_Schiller") & (x["origin"] == "normal"),
    "origin",
] = "normal_adjacent"
# Normal samples from Laughney_Massague should be normal adjacent
adata.obs.loc[
    lambda x: x["dataset"].str.contains("Laughney_Massague")
    & (x["origin"] == "normal"),
    "origin",
] = "normal_adjacent"

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
# ### Add "study" in addition to dataset

# %%
adata.obs["study"] = ["_".join(x.split("_")[:3]) for x in adata.obs["dataset"]]

# %%
adata.obs["study"].unique()

# %% [markdown]
# # Add sequencing platforms

# %%
platform_metadata = pd.read_csv(path_platform_metadata)

# %%
platform_metadata_by_cell = (
    adata.obs.loc[:, ["dataset"]]
    .reset_index(drop=False)
    .merge(platform_metadata, on="dataset", how="left")
    .set_index("index")
)

# %%
adata.obs["platform"] = platform_metadata_by_cell["platform"]
adata.obs["platform_fine"] = platform_metadata_by_cell["platform_fine"]

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

# %% [markdown]
# ### cell_type_coarse

# %%
cell_type_coarse_map = {
    "Alveolar cell type 1": "Epithelial cell",
    "Alveolar cell type 2": "Epithelial cell",
    "B cell": "B cell",
    "B cell dividing": "B cell",
    "Ciliated": "Epithelial cell",
    "Club": "Epithelial cell",
    "DC mature": "cDC",
    "Endothelial cell arterial": "Endothelial cell",
    "Endothelial cell capillary": "Endothelial cell",
    "Endothelial cell lymphatic": "Endothelial cell",
    "Endothelial cell venous": "Endothelial cell",
    "Fibroblast adventitial": "Stromal",
    "Fibroblast alveolar": "Stromal",
    "Fibroblast peribronchial": "Stromal",
    "Macrophage": "Macrophage/Monocyte",
    "Macrophage alveolar": "Macrophage/Monocyte",
    "Mast cell": "Mast cell",
    "Mesothelial": "Stromal",
    "Monocyte classical": "Macrophage/Monocyte",
    "Monocyte non-classical": "Macrophage/Monocyte",
    "NK cell": "NK cell",
    "Neutrophils": "Neutrophils",
    "Pericyte": "Stromal",
    "Plasma cell": "Plasma cell",
    "Plasma cell dividing": "Plasma cell",
    "ROS1+ healthy epithelial": "Epithelial cell",
    "Smooth muscle cell": "Stromal",
    "T cell CD4": "T cell",
    "T cell CD8": "T cell",
    "T cell dividing": "T cell",
    "T cell regulatory": "T cell",
    "Tumor cells": "Epithelial cell",
    "cDC1": "cDC",
    "cDC2": "cDC",
    "myeloid dividing": "Macrophage/Monocyte",
    "pDC": "pDC",
    "stromal dividing": "Stromal",
    "transitional club/AT2": "Epithelial cell",
}

# %%
adata.obs["cell_type_coarse"] = [
    cell_type_coarse_map[x] for x in adata.obs["cell_type"]
]

# %% [markdown]
# ### cell_type_major
#
# an annotation that is more fine-grained than coarse, more coarse grained than fine, and ignore some special categories, such as dividing cells. 
# This may be useful for some downstream analyses. 

# %%
cell_type_major_map = {
    "Alveolar cell type 1": "Alveolar cell type 1",
    "Alveolar cell type 2": "Alveolar cell type 2",
    "B cell": "B cell",
    "B cell dividing": "other",
    "Ciliated": "Ciliated",
    "Club": "Club",
    "DC mature": "DC mature",
    "Endothelial cell arterial": "Endothelial cell",
    "Endothelial cell capillary": "Endothelial cell",
    "Endothelial cell lymphatic": "Endothelial cell",
    "Endothelial cell venous": "Endothelial cell",
    "Fibroblast adventitial": "Stromal",
    "Fibroblast alveolar": "Stromal",
    "Fibroblast peribronchial": "Stromal",
    "Macrophage": "Macrophage",
    "Macrophage alveolar": "Macrophage alveolar",
    "Mast cell": "Mast cell",
    "Mesothelial": "Stromal",
    "Monocyte classical": "Monocyte",
    "Monocyte non-classical": "Monocyte",
    "NK cell": "NK cell",
    "Neutrophils": "Neutrophils",
    "Pericyte": "Stromal",
    "Plasma cell": "Plasma cell",
    "Plasma cell dividing": "other",
    "ROS1+ healthy epithelial": "other",
    "Smooth muscle cell": "Stromal",
    "T cell CD4": "T cell CD4",
    "T cell CD8": "T cell CD8",
    "T cell dividing": "other",
    "T cell regulatory": "T cell regulatory",
    "Tumor cells": "Tumor cells",
    "cDC1": "cDC1",
    "cDC2": "cDC2",
    "myeloid dividing": "other",
    "pDC": "pDC",
    "stromal dividing": "other",
    "transitional club/AT2": "transitional club/AT2",
}

# %%
adata.obs["cell_type_major"] = [cell_type_major_map[x] for x in adata.obs["cell_type"]]

# %% [markdown]
# # Dataset overview

# %% [markdown]
# ### Overview of celltypes

# %%
sc.pl.umap(adata, color="total_counts", vmax=10000)

# %%
sc.set_figure_params(figsize=(10, 10))
sc.pl.umap(
    adata,
    color=["cell_type_coarse", "cell_type", "cell_type_major", "platform"],
    legend_loc="on data",
    legend_fontoutline=2,
    legend_fontsize=8,
    size=1,
)

# %%
sc.set_figure_params(figsize=(6, 6))
sc.pl.umap(
    adata,
    color=[
        "condition",
        "origin",
        "dataset",
        "sex",
        "ever_smoker",
        "uicc_stage",
        "tumor_stage",
    ],
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

# %%

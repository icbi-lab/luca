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
import scvi

sc.settings.set_figure_params(figsize=(4, 4))
import pandas as pd
from nxfvars import nxfvars
from scanpy_helpers.integration import (
    sanitize_adata,
    validate_adata,
    aggregate_duplicate_gene_symbols,
)
import scipy.sparse
import anndata
import sklearn
from threadpoolctl import threadpool_limits
import matplotlib.pyplot as plt
import scanpy_helpers as sh
from scanpy_helpers.annotation import AnnotationHelper
import numpy as np
from tqdm import tqdm
from pathlib import Path

# %% [markdown]
# ## Get Input data

# %%
threadpool_limits(nxfvars.get("cpus", 16))
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")
dataset_path = nxfvars.get(
    "dataset_path",
    "../../data/30_downstream_analyses/02_integrate_into_atlas/",
)
reference_scanvi = sc.read_h5ad(
    nxfvars.get(
        "reference_scanvi_h5ad",
        "../../data/20_build_atlas/annotate_datasets/35_final_atlas/full_atlas_hvg_integrated_scvi_integrated_scanvi.h5ad",
    )
)
reference_atlas = sc.read_h5ad(
    nxfvars.get(
        "reference_atlas",
        "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
    )
)

# %%
dataset_paths = list(Path(dataset_path).glob("**/*_integrated.h5ad"))
datasets = {
    p.stem.replace("_integrated", ""): sc.read_h5ad(p) for p in tqdm(dataset_paths)
}

# %% [markdown]
# ## Fix Metadata

# %%
datasets["UKIM-V-2"].obs["dataset"] = "UKIM-V-2"
datasets["UKIM-V-2"].obs["patient"] = [
    f"UKIM-V-2_{p}" for p in datasets["UKIM-V-2"].obs["patient"]
]
datasets["UKIM-V-2"].obs["tissue"] = "lung"
datasets["UKIM-V-2"].obs["tumor_stage"] = np.nan
datasets["UKIM-V-2"].obs.loc[
    datasets["UKIM-V-2"].obs["uicc_stage"].isin(["I", "II", "IA"]), "tumor_stage"
] = "early"
datasets["UKIM-V-2"].obs["sample"] = "UKIM-V-2_" + datasets["UKIM-V-2"].obs[
    "sample"
].astype(str)
datasets["UKIM-V-2"].obs.loc[
    datasets["UKIM-V-2"]
    .obs["uicc_stage"]
    .isin(["III", "IV", "IIIA", "IIIB", "III or IV"]),
    "tumor_stage",
] = "late"
assert np.all(
    pd.isnull(datasets["UKIM-V-2"].obs["uicc_stage"])
    == pd.isnull(datasets["UKIM-V-2"].obs["tumor_stage"])
)

# %% [markdown] tags=[]
# ## Merge query and reference

# %%
for dataset_id, adata in datasets.items():
    adata.obs["dataset"] = dataset_id
    validate_adata(adata)

# %%
for dataset in datasets.values():
    assert dataset.shape[1] == datasets["UKIM-V-2"].shape[1], "shape mismatch"

# %%
query_merged = anndata.concat(
    datasets.values(),
    index_unique="-",
)

# %%
reference_atlas.obsm["X_scANVI"] = reference_scanvi.obsm["X_scANVI"]

# %%
adata_full = anndata.concat(
    [reference_atlas, query_merged],
    merge="first",
    uns_merge="first",
    join="outer",
)
assert adata_full.var_names.is_unique
assert adata_full.obs_names.is_unique

# %%
# raw does not get concatenated... let's recompute it
adata_full_raw = adata_full.copy()
sc.pp.normalize_total(adata_full_raw)
sc.pp.log1p(adata_full_raw)
adata_full.raw = adata_full_raw

# %%
# Excluding Maier_Merad_2020 which is a subset of Leader_Merad_2021
adata_full = adata_full[adata_full.obs["dataset"] != "Maier_Merad_2020", :]
adata_full = adata_full[adata_full.obs["doublet_status"] != "doublet", :].copy()

# %%
sc.pp.neighbors(adata_full, use_rep="X_scANVI")


# %%
# initalize umap with the original coordinates. Missing values (=query) are initialized with random values.
init_pos_df = pd.DataFrame(
    reference_scanvi.obsm["X_umap"], index=reference_scanvi.obs_names
).reindex(adata_full.obs_names)
for col in init_pos_df.columns:
    na_mask = init_pos_df[col].isnull()
    init_pos_df.loc[na_mask, col] = np.random.uniform(
        np.min(init_pos_df[col]), np.max(init_pos_df[col]), size=np.sum(na_mask)
    )

# %%
sc.tl.umap(adata_full, init_pos=init_pos_df.values)

# %% [markdown]
# ## Automated cell-type annotation
#
# The annotation by scANVI doesn't look very good. I, therefore, implemented an
# approach based on nearest neighbor majority voting.

# %%
col = "cell_type_tumor"
# I want the prediction in a separate column...
sh.annotation.classify_cell_types_nearest_neighbors(
    adata_full,
    col,
    mask_reference=~adata_full.obs["dataset"].isin(datasets.keys()),
    mask_query=adata_full.obs["dataset"].isin(datasets.keys()),
    key_added=f"_{col}_predicted",
)
# ... but also merged back into the original column
adata_full.obs.loc[~adata_full.obs[f"_{col}_predicted"].isnull(), col] = adata_full.obs[
    f"_{col}_predicted"
]

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(adata_full, color="dataset", groups=list(datasets.keys()), size=4)

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(adata_full, color=[col, f"_{col}_predicted"], size=2)

# %% [markdown]
# ### Re-annotate more coarse-grained categories

# %%
ah = sh.annotation.AnnotationHelper()

# %%
adata_full.obs["cell_type"] = adata_full.obs["cell_type_tumor"].astype(str)
adata_full.obs.loc[
    adata_full.obs["cell_type_tumor"].str.startswith("Tumor"), "cell_type"
] = "Tumor cells"
adata_full.obs["cell_type_tumor"] = adata_full.obs["cell_type_tumor"].str.replace("LSCC", "LUSC")

# %%
with plt.rc_context({"figure.figsize": (5, 5)}):
    ah.annotate_cell_types(
        adata_full,
        {
            "Alveolar cell type 1": ["Alveolar cell type 1"],
            "Alveolar cell type 2": ["Alveolar cell type 2"],
            "B cell": ["B cell"],
            "Ciliated": ["Ciliated"],
            "Club": ["Club"],
            "DC mature": ["DC mature"],
            "Endothelial cell": ["Endothelial cell", "Endothelial cell lymphatic"],
            "Goblet": ["Goblet"],
            "Macrophage": ["Macrophage"],
            "Macrophage FABP4+": ["Macrophage FABP4+"],
            "Mast cell": ["Mast cell"],
            "Monocyte": ["Monocyte"],
            "NK cell": ["NK cell"],
            "Neutrophils": ["Granulocytes"],
            "Plasma cell": ["Plasma cell"],
            "Stromal": [
                "Fibroblast",
                "Fibroblast adventitial",
                "Fibroblast alveolar",
                "Smooth muscle cell",
                "Mesothelial",
                "Pericyte",
            ],
            "T cell CD4": ["T cell CD4"],
            "T cell CD8": ["T cell CD8"],
            "T cell regulatory": ["T cell regulatory"],
            "Tumor cells": ["Tumor cells"],
            "cDC1": ["cDC1"],
            "cDC2": ["cDC2"],
            "pDC": ["pDC"],
        },
        key_added="cell_type_major",
        column="cell_type",
    )

# %%
with plt.rc_context({"figure.figsize": (5, 5)}):
    ah.annotate_cell_types(
        adata_full,
        {
            "Epithelial cell": [
                "Alveolar cell type 1",
                "Alveolar cell type 2",
                "Ciliated",
                "Club",
                "Tumor cells",
                "Goblet",
                "ROS1+ healthy epithelial",
            ],
            "B cell": ["B cell", "B cell dividing"],
            "cDC": ["DC mature", "cDC1", "cDC2"],
            "Endothelial cell": ["Endothelial cell", "Endothelial cell lymphatic"],
            "Macrophage/Monocyte": [
                "Macrophage",
                "Macrophage FABP4+",
                "myeloid dividing",
                "Monocyte",
            ],
            "Mast cell": ["Mast cell"],
            "NK cell": ["NK cell"],
            "Neutrophils": ["Granulocytes"],
            "Plasma cell": ["Plasma cell"],
            "Stromal": [
                "Fibroblast",
                "Fibroblast adventitial",
                "Fibroblast alveolar",
                "Smooth muscle cell",
                "Mesothelial",
                "Pericyte",
            ],
            "T cell": [
                "T cell CD4",
                "T cell CD8",
                "T cell regulatory",
                "T cell dividing",
                "other (T assoc.)",
            ],
            "pDC": ["pDC"],
        },
        key_added="cell_type_coarse",
        column="cell_type",
    )

# %% [markdown]
# ## Update dataset-level metadata

# %%
adata_full.obs["study"] = ["_".join(x.split("_")[:3]).replace("-2", "") for x in adata_full.obs["dataset"]]

# %%
adata_full.obs["study"].unique()

# %% [markdown]
# ## Write output
#

# %%
adata_full.write_h5ad(f"{artifact_dir}/full_atlas_merged.h5ad")

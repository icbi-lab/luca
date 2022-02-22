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
import pandas as pd
from nxfvars import nxfvars
import numpy as np
from scanpy_helpers import compare_groups

sc.settings.set_figure_params(figsize=(5, 5))

# %% [markdown]
# # Load input data

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")
adata_in = nxfvars.get(
    "adata_in",
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
)
adata_in_cpdb = nxfvars.get(
    "adata_in_cpdb",
    "../../data/30_downstream_analyses/cell2cell/cpdb_h5ad/artifacts/adata_cpdb.h5ad",
)
stratification_csv = nxfvars.get(
    "stratification_csv",
    "../../data/30_downstream_analyses/stratify_patients/artifacts/patient_stratification.csv",
)
# run one of the comparisons defined in the list below
comparison = nxfvars.get("comparison", "tumor_normal")
cpus = nxfvars.get("cpus", 32)

# %%
patient_strat = pd.read_csv(stratification_csv, index_col=0)

# %%
adata = sc.read_h5ad(adata_in)

# %%
adata_cpdb = sc.read_h5ad(adata_in_cpdb)

# %%
sc.pl.umap(adata, color=["cell_type_major"], wspace=1)

# %% [markdown] tags=[]
# # Prepare datasets
# ## adata tumor/normal
#  * includes tumor primary and normal/normal adjacent samples
#  * normal and normal adjacent are merged into "normal"

# %%
adata_tumor_normal = adata[
    adata.obs["origin"].isin(["normal_adjacent", "normal", "tumor_primary"]),
    :,
].copy()
adata_tumor_normal.obs["origin"] = [
    "normal" if "normal" in x else x for x in adata_tumor_normal.obs["origin"]
]
adata_tumor_normal.obs["origin"].unique()

# %% [markdown]
# ## adata primary tumor
#  * includes only primary tumor samples
#  * excludes datasets which only include a single cell-types

# %%
adata_primary_tumor = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"]),
    :,
]

# %%
adata_primary_tumor.obs = (
    adata_primary_tumor.obs.reset_index()
    .merge(
        patient_strat,
        how="left",
    )
    .set_index("index")
)

# %%
# Create additional variables based on others
adata_primary_tumor.obs["infiltration_status"] = [
    {"-": "immune_low", "T": "immune_high", "B": "immune_high", "M": "immune_high"}.get(
        i, np.nan
    )
    for i in adata_primary_tumor.obs["immune_infiltration"]
]
adata_primary_tumor.obs["infiltration_type"] = [
    np.nan if i == "-" else i for i in adata_primary_tumor.obs["immune_infiltration"]
]

# %%
# LSCCness defined as the fraction of tumor cells of all tumor cells that are assigned to one of the LSCC clusters.
lsccness = (
    adata_primary_tumor.obs.groupby("patient")
    .apply(
        lambda x: np.sum(x["cell_type_tumor"].str.contains("LSCC"))
        / np.sum(x["cell_type_major"] == "Tumor cells")
    )
    .reset_index(name="lsccness")
    .merge(adata_primary_tumor.obs.loc[:, "patient"].reset_index())
    .set_index("index")
)

# %%
adata_primary_tumor.obs["lsccness"] = lsccness["lsccness"]

# %%
adata_cpdb.obs["lsccness"] = (
    adata_primary_tumor.obs.loc[:, ["sample", "lsccness"]]
    .drop_duplicates()
    .assign(idx=lambda x: x["sample"].str.lower())
    .set_index("idx")["lsccness"]
)

# %% [markdown]
# # List of comparisons

# %%
comparisons = {
    "tumor_normal": {
        "dataset": adata_tumor_normal,
        "dataset_cpdb": adata_cpdb,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "origin",
        "lm_covariate_str": "+ patient",  # paired analysis by patient
        "contrasts": "Treatment('normal')",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
    "infiltration_status": {
        "dataset": adata_primary_tumor[
            adata_primary_tumor.obs["infiltration_status"].isin(
                ["immune_low", "immune_high"]
            ),
            :,
        ],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "infiltration_status",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Treatment('immune_low')",
    },
    "infiltration_type": {
        "dataset": adata_primary_tumor[
            adata_primary_tumor.obs["infiltration_type"].isin(["T", "B", "M"]), :
        ],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "infiltration_type",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Sum",
    },
    "patient_immune_infiltration": {
        "dataset": adata_primary_tumor,
        "dataset_cpdb": adata_cpdb,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "immune_infiltration",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Sum",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
    "patient_immune_infiltration_condition": {
        "dataset": adata_primary_tumor[
            ~adata_primary_tumor.obs["lsccness"].isnull(), :
        ],
        "dataset_cpdb": adata_cpdb[~adata_cpdb.obs["lsccness"].isnull(), :],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient", "lsccness"],
        "column_to_test": "immune_infiltration",
        "lm_covariate_str": "+ dataset + lsccness",
        "contrasts": "Sum",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
    "patient_immune_infiltration_treatment_coding": {
        "dataset": adata_primary_tumor,
        "dataset_cpdb": adata_cpdb,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "immune_infiltration",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Treatment('-')",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
    "patient_immune_infiltration_treatment_coding_condition": {
        "dataset": adata_primary_tumor[
            ~adata_primary_tumor.obs["lsccness"].isnull(), :
        ],
        "dataset_cpdb": adata_cpdb[~adata_cpdb.obs["lsccness"].isnull(), :],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient", "lsccness"],
        "column_to_test": "immune_infiltration",
        "lm_covariate_str": "+ dataset + lsccness",
        "contrasts": "Treatment('-')",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
    "patient_immune_infiltration_treatment_coding_condition2": {
        "dataset": adata_primary_tumor[
            adata_primary_tumor.obs["condition"].isin(["LUAD", "LSCC"]), :
        ],
        "dataset_cpdb": adata_cpdb,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient", "condition"],
        "column_to_test": "immune_infiltration",
        "lm_covariate_str": "+ dataset + condition",
        "contrasts": "Treatment('-')",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
    "luad_lscc": {
        "dataset": adata_primary_tumor[
            adata_primary_tumor.obs["condition"].isin(["LUAD", "LSCC"]), :
        ],
        "dataset_cpdb": adata_cpdb,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient", "tumor_stage"],
        "column_to_test": "condition",
        "lm_covariate_str": "+ dataset + tumor_stage",
        "contrasts": "Treatment('LUAD')",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
    "early_advanced": {
        "dataset": adata_primary_tumor,
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient"],
        "column_to_test": "tumor_stage",
        "lm_covariate_str": "+ dataset",
        "contrasts": "Treatment('early')",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
    "early_advanced_condition": {
        "dataset": adata_primary_tumor[
            ~adata_primary_tumor.obs["lsccness"].isnull(), :
        ],
        "cell_type_column": "cell_type_major",
        "pseudobulk_group_by": ["dataset", "patient", "lsccness"],
        "column_to_test": "tumor_stage",
        "lm_covariate_str": "+ dataset + lsccness",
        "contrasts": "Treatment('early')",
        "tools": ["dorothea", "progeny", "cytosig"],
    },
}

# %% [markdown]
# # Run differential signature analysis

# %%
comparison_config = comparisons[comparison]

# %%
datasets = compare_groups.prepare_dataset(comparison, n_jobs=cpus, **comparison_config)

# %%
results = compare_groups.compare_signatures(
    comparison, datasets, n_jobs=cpus, **comparison_config
)

# %%
if "dataset_cpdb" in comparison_config:
    results["cpdb"] = compare_groups.compare_cpdb(
        comparison, n_jobs=cpus, **comparison_config
    )

# %%
for tool, results in results.items():
    results.to_csv(f"{artifact_dir}/{comparison}_{tool}.tsv", sep="\t")

# %%

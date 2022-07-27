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
import numpy.testing as npt

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
    "../../data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
)
path_neutrophil_atlas = nxfvars.get(
    "neutrophil_atlas",
    "../../data/30_downstream_analyses/neutrophils/subclustering/artifacts/full_atlas_neutrophil_clusters.h5ad",
)

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
# ## Matrix layers
#
#  * all genes
#  * sparse
#  * `.raw` -> UMI counts or read counts, respectively
#  * "final" -> suitable for visualization in cellxgene (log-norm, lenght-scaled for SS2) 
#  * `layers["counts_length_scaled"]` for UMI counts and length-scaled SS2 counts, respectively. 
#

# %%
extended_atlas

# %%
new_raw = sc.AnnData(core_atlas.layers["raw_counts"], var=core_atlas.var)
core_atlas.layers["counts_length_scaled"] = core_atlas.X
core_atlas.X = core_atlas.raw.X
core_atlas.raw = new_raw
del core_atlas.layers["raw_counts"]

# %%
# small consistency check.
for atlas in [core_atlas, extended_atlas]:
    npt.assert_array_equal(
        atlas[
            ["SRR10777354-13"], ["A1BG", "A2M", "AAK1", "AAMP", "ABCA1", "ABCD1"]
        ].X.A,
        [[20.0, 18.0, 20.0, 20.0, 93.0, 6.0]],
    )
    npt.assert_array_equal(
        atlas[["SRR10777354-13"], ["A1BG", "A2M", "AAK1", "AAMP", "ABCA1", "ABCD1"]]
        .layers["raw_counts"]
        .A,
        [[3.0, 66.0, 142.0, 12.0, 522.0, 15.0]],
    )

# %% [markdown]
# ## Reannotate with ontologies

# %% [markdown]
# ### Sequencing platform / assay_ontology_term_id
#
# If there is not an exact match for the assay, clarifying text MAY be enclosed in parentheses and appended to the most accurate term. For example, the sci-plex assay could be curated as "EFO:0010183 (sci-plex)".

# %%
extended_atlas.obs["platform_fine"].value_counts()

# %%
platform_map = {
    "10x_3p_v2": "EFO:0009899",
    "10x_3p_v3": "EFO:0009922",
    "10x_5p_v1": "EFO:0011025",
    "Smart-seq2": "EFO:0008931",
    "DropSeq": "EFO:0008722",
    "Singleron": "EFO:0030031",
    "BD-Rhapsody": "EFO:0010183 (BD Rhapsody)",  # single-cell library construction. Cannot find BD Rhapsody in ontology
    "InDrop": "EFO:0008780",
}

# %% [markdown]
# ### cell types / cell_type_ontology_term_id

# %%
# based on the fines cell-type annotation that does not use specific subclusters
cell_type_map = {
    "Alveolar cell type 1": "CL:0002062",
    "Alveolar cell type 2": "CL:0002063",
    "B cell": "CL:0000236",
    "B cell dividing": "CL:0000236",
    "Ciliated": "CL:0005012",
    "Club": "CL:0000158",
    "DC mature": "CL:0000451",
    "Endothelial cell arterial": "CL:1001568",
    "Endothelial cell capillary": "CL:0002144",
    "Endothelial cell lymphatic": "CL:0002138",
    "Endothelial cell venous": "CL:0002543",
    "Fibroblast adventitial": "CL:0002553",
    "Fibroblast alveolar": "CL:0002553",
    "Fibroblast peribronchial": "CL:2000093",
    "Macrophage": "CL:0000235",
    "Macrophage alveolar": "CL:0000583",
    "Mast cell": "CL:0000097",
    "Mesothelial": "CL:0000077",
    "Monocyte classical": "CL:0000860",
    "Monocyte non-classical": "CL:0000875",
    "NK cell": "CL:0000623",
    "NK cell dividing": "CL:0000623",
    "Neutrophils": "CL:0000775",
    "Pericyte": "CL:0000669",
    "Plasma cell": "CL:0000786",
    "Plasma cell dividing": "CL:0000786",
    "ROS1+ healthy epithelial": "CL:0000082",
    "Smooth muscle cell": "CL:0000192",
    "T cell CD4": "CL:0000624",
    "T cell CD4 dividing": "CL:0000624",
    "T cell CD8 activated": "CL:0000625",
    "T cell CD8 dividing": "CL:0000625",
    "T cell CD8 effector memory": "CL:0000625",
    "T cell CD8 naive": "CL:0000625",
    "T cell CD8 terminally exhausted": "CL:0000625",
    "T cell NK-like": "CL:0000625",
    "T cell regulatory": "CL:0000815",
    "Tumor cells": "CL:0001064",
    "cDC1": "CL:0000990",
    "cDC2": "CL:0002399",
    "myeloid dividing": "CL:0000763",
    "pDC": "CL:0000784",
    "stromal dividing": "CL:0000499",
    "transitional club/AT2": "CL:0000082",
}
assert "" not in cell_type_map.values()


# %% [markdown]
# ### age / development_stage_ontology_term_id
#
# If unavailable, must be `unknown`

# %%
def get_stage(age):
    age = int(age)
    if age < 20:
        raise NotImplementedError()
    elif age <= 79:
        return f"HsapDv_0000{age + 94}"
    elif age <= 99:
        return f"HsapDv_0000{age + 126}"
    else:
        raise NotImplementedError()


# %% [markdown]
# ### condition / disease_ontology_term_id
#
#  This MUST be a MONDO term or "PATO:0000461" for normal or healthy. 

# %%
extended_atlas.obs["condition"].value_counts()

# %%
condition_map = {
    "LUAD": "MONDO:0005061",
    "non-cancer": "PATO:0000461",
    "LUSC": "MONDO:0005097",
    "NSCLC NOS": "MONDO:0005233",
    "COPD": "MONDO:0005002",
}


# %% [markdown]
# ### ethnicity_ontology_term_id
#
# for Homo sapiens, this MUST be either a HANCESTRO term or "unknown" if unavailable. 

# %% [markdown]
# ### is_primary_data

# %%
def is_primary_data(dataset):
    return "UKIM" in dataset


# %% [markdown]
# ### organism_ontology_term_id
#
# `NCBITaxon:9606` for human

# %% [markdown]
# ### sex_ontology_term_id
# This MUST be a child of PATO:0001894 for phenotypic sex or "unknown" if unavailable.

# %%
extended_atlas.obs["sex"].value_counts()

# %%
sex_map = {
    "male": "PATO_0000384",
    "female": "PATO_0000383",
}

# %% [markdown]
# ### tissue_ontology_term_id

# %%
extended_atlas.obs["tissue"].value_counts()

# %%
tissue_map = {
    "lung": "UBERON:0002048",
    "lymph_node": "UBERON:0000029",
    "pleura/effusion": "UBERON:0000175",
    "brain": "UBERON:0000955",
    "liver": "UBERON:0002107",
    "adrenal": "UBERON:0018303",
}

# %% [markdown]
# ## Gene annotations
#  * GENCODE v38/Ensmbl 104
#  
# A gene symbol to ENSEMBL map was generated based on the gencode.v38 GTF file (see commented code below).
# If an ensembl id was missing for a gene symbol, we manually complemented it with information from genecards.org. Duplicate entries were manually resolved using genecards.org. 

# %%
# gtf = gtfparse.read_gtf(
#     "../../data/10_references/gencode.v38.primary_assembly.annotation.gtf.gz",
#     usecols=["gene_id", "gene_name"],
# ).drop_duplicates()
# var2 = extended_atlas.var.join(gtf.set_index("gene_name"))
# var2.loc[:, ["gene_id"]].reset_index().rename(columns={"index": "gene_symbol"}).to_csv("../../tables/symbol_to_ensembl.csv")

# %%
gene2ensembl = pd.read_csv("../../tables/symbol_to_ensembl.csv", index_col=0)
gene2ensembl["gene_id"] = (
    gene2ensembl["gene_id"].str.strip().str.replace("\.\d+$", "", regex=True)
)
gene2ensembl["gene_id2"] = [
    gene_id if not pd.isnull(gene_id) else gene_symbol
    for gene_id, gene_symbol in zip(
        gene2ensembl["gene_id"], gene2ensembl["gene_symbol"]
    )
]
gene2ensembl = gene2ensembl.loc[lambda x: ~x["gene_id2"].str.contains("_PAR_Y"), :]

# %%
extended_atlas.var = extended_atlas.var.join(
    gene2ensembl.set_index("gene_symbol"), how="left"
)

# %%
extended_atlas.var.set_index("gene_id2", inplace=True, drop=False)

# %%
pd.set_option("display.max_rows", 100)
extended_atlas.var.loc[lambda x: x["gene_id"].isnull()]

# %% [markdown]
# ## Gene metadata / var and raw.var

# %%
extended_atlas.var

# %% [markdown]
# ### feature_biotype
#
# must be `gene` unless it's a spike-in

# %%
extended_atlas.var["feature_biotype"] = "gene"

# %% [markdown]
# ### Index
#
# must be set to ENSEMBL gene names

# %% [markdown]
# ### feature_is_filtered
# * Curators MUST annotate the following column only in the var dataframe. This column MUST NOT be present in raw.var
# * This MUST be True if the feature was filtered out in the final matrix (X) but is present in the raw matrix

# %%
extended_atlas.var["feature_is_filtered"] = False

# %% [markdown]
# ## Metadata / uns

# %%
extended_atlas.uns["schema_version"] = "2.0.0"
extended_atlas.uns[
    "title"
] = "The single-cell lung cancer atlas (LuCA) -- extended atlas"
extended_atlas.uns["X_normalization"] = (
    "Log-normalized length-scaled counts (Smart-seq2) or UMI counts (other protocols). "
    "Normalization was performed using sc.pp.normalize_total() followed by sc.pp.log1p()"
)
extended_atlas.uns["batch_condition"] = ["dataset", "patient", "sample"]
extended_atlas.uns["default_embedding"] = "X_umap"

# %% [markdown]
# ## Export h5ads

# %% [markdown]
# ## Validate h5ads

# %%
# TODO
# !cellxgene-schema validate {path_core_atlas}

# %%

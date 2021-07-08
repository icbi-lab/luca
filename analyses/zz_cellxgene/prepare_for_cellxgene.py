# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python [conda env:.conda-scanpy_2020-12]
#     language: python
#     name: conda-env-.conda-scanpy_2020-12-py
# ---

# %% [markdown]
# ### Prepare datasets for cellxgene
#
# -> all genes in `X`, normalized values in `X`. 

# %% [markdown]
# ## Maynard

# %%
import scanpy as sc

# %%
maynard_raw = sc.read_h5ad(
    "../../data/10_public_datasets/Maynard_Bivona_2020_NSCLC/h5ad_raw/maynard2020.h5ad"
)

# %%
maynard = sc.read_h5ad("../../data/30_annotate_scrnaseq_data/maynard_annotated.h5ad")

# %%
maynard_raw_subset = maynard_raw[maynard.obs_names, :]

# %%
maynard_cellxgene = sc.AnnData(
    var=maynard_raw_subset.var,
    X=maynard_raw_subset.layers["tpm"],
    obs=maynard.obs,
    obsm=maynard.obsm,
)

# %%
sc.pp.log1p(maynard_cellxgene)

# %%
sc.pl.umap(maynard_cellxgene, color="STING1")

# %%
# !mkdir -p ../../data/zz_cellxgene
maynard_cellxgene.write_h5ad("../../data/zz_cellxgene/maynard2020.h5ad", compression="lzf")

# %% [markdown]
# ## Lambrechts

# %%
lambrechts_raw = sc.read_h5ad(
    "../../data/10_public_datasets/Lambrechts_2018_LUAD/E-MTAB-6653/h5ad_raw/lambrechts_2018_luad_6653.h5ad"
)

# %%
lambrechts = sc.read_h5ad(
    "../../data/30_annotate_scrnaseq_data/lambrechts_annotated.h5ad"
)

# %%
lambrechts_raw_subset = lambrechts_raw[lambrechts.obs_names, :]

# %%
lambrechts_cellxgene = sc.AnnData(
    var=lambrechts_raw_subset.var,
    X=lambrechts_raw_subset.X,
    obs=lambrechts.obs,
    obsm=lambrechts.obsm,
)

# %%
sc.pp.normalize_total(lambrechts_cellxgene)
sc.pp.log1p(lambrechts_cellxgene)

# %%
lambrechts_cellxgene.write_h5ad("../../data/zz_cellxgene/lambrechts2020_6653.h5ad", compression="lzf")

# %% [markdown]
# ## Integrated all

# %%
integrated_all = sc.read_h5ad("../../data/50_integrate_scrnaseq_data/integrated_merged_all.h5ad")

# %%
integrated_all_cellxgene = sc.AnnData(
    var=integrated_all.raw.var,
    X=integrated_all.raw.X,
    obs=integrated_all.obs,
    obsm=integrated_all.obsm
)

# %%
integrated_all_cellxgene.write_h5ad("../../data/zz_cellxgene/scanvi_integrated_all.h5ad", compression="lzf")

# %% [markdown]
# ## Integrated tumor samples

# %%
integrated_cancer = sc.read_h5ad("../../data/50_integrate_scrnaseq_data/integrated_merged_nsclc_heterogeneity.h5ad")

# %%
integrated_cancer_cellxgene = sc.AnnData(
    var=integrated_cancer.raw.var,
    X=integrated_cancer.raw.X,
    obs=integrated_cancer.obs,
    obsm=integrated_cancer.obsm
)

# %%
integrated_cancer_cellxgene.write_h5ad("../../data/zz_cellxgene/scanvi_integrated_tumor_samples.h5ad", compression="lzf")

# %%

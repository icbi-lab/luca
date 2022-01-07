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
from scanpy_helpers.annotation import AnnotationHelper
import scanpy_helpers as sh
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf


# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
adata = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad"
)

# %%
adata_t = adata[adata.obs["cell_type_coarse"] == "Tumor cells", :].copy()
adata_n = adata[adata.obs["cell_type_coarse"] == "Neutrophils", :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_n, use_rep="X_scANVI")

# %%
sc.pl.umap(
    adata_n,
    color=["cell_type", "origin", "condition", "tumor_stage", "sex", "dataset"],
    wspace=0.5,
)

# %% [markdown]
# ### Compare cytosig and progeny

# %%
tumor_normal_progeny = pd.read_csv(
    "/home/sturm/Downloads/differential_signature_origin_progeny.tsv",
    sep="\t",
    index_col=0,
)
tumor_normal_cytosig = pd.read_csv(
    "/home/sturm/Downloads/differential_signature_origin_cytosig.tsv",
    sep="\t",
    index_col=0,
)

# %%
tumor_normal_progeny.loc[lambda x: x["cell_type"] == "Neutrophils", :].query(
    "pvalue < 0.1"
)

# %%
tumor_normal_cytosig.loc[lambda x: x["cell_type"] == "Neutrophils", :].query(
    "pvalue < 0.1"
)

# %% [markdown]
# ### Compare DE genes using pair plots

# %%
top_up_normal = [
    "S100A8",
    "IFITM2",
    "ACTB",
    "RGS2",
    "SELL",
    "S100A9",
    "C5AR1",
    "MNDA",
    "PTGS2",
    "ARHGDIB",
    "H3-3A",
    "FPR1",
    "CXCR2",
    "CSF3R",
    "CNN2",
    "SRGN",
    "CTSS",
    "TPM4",
    "TAGLN2",
    "VNN2",
]

# %%
top_up_tumor = [
    # "CCL3L1",
    "C15orf48",
    "CXCR4",
    "VEGFA",
    "CD83",
    "TPI1",
    "CCL4L2",
    "CCRL2",
    "IER3",
    "BRI3",
    "SQSTM1",
    "JUN",
    "PHACTR1",
    "GAPDH",
    "CSTB",
    "LGALS3",
    "ENO1",
    "NFKBIA",
    "FCER1G",
    "WSB1",
]

# %%
tumor_vs_normal = ["PTGS2", "SELL", "CXCR2", "VEGFA", "OLR1", "CXCR4"]

# %%
tmp_adata = adata_n[
    adata_n.obs["origin"].isin(["normal_adjacent", "normal", "tumor_primary"]),
    :,
]
tmp_adata.obs["origin"] = [
    "normal" if "normal" in x else x for x in tmp_adata.obs["origin"]
]
adata_tumor_normal = sh.pseudobulk.pseudobulk(
    tmp_adata, groupby=["dataset", "patient", "origin"]
)
sc.pp.normalize_total(adata_tumor_normal, target_sum=1000)
sc.pp.log1p(adata_tumor_normal)

# %%
adata_tumor_normal.obs.query("dataset == 'UKIM-V'")

# %%
sc.pl.matrixplot(
    adata_tumor_normal, var_names=top_up_normal, groupby="origin", title="top up normal"
)

# %%
sc.pl.matrixplot(
    adata_tumor_normal, var_names=top_up_tumor, groupby="origin", title="top up tumor"
)

# %%
sh.pairwise.plot_paired(
    adata_tumor_normal, groupby="origin", paired_by="patient", var_names=tumor_vs_normal
)

# %%
sh.pairwise.plot_paired(
    adata_tumor_normal, groupby="origin", var_names=tumor_vs_normal, hue="dataset", show_legend=False
)

# %% [markdown]
# ## neutro fractions by tumor type

# %%
ct_fractions = (
    adata[adata.obs["origin"] == "tumor_primary", :]
    .obs.groupby(["dataset", "patient", "condition"])["cell_type"]
    .value_counts(normalize=True)
    .reset_index()
)

# %%
ct_fractions = ct_fractions.rename(
    columns={"level_3": "cell_type", "cell_type": "fraction"}
)

# %%
neutro_fractions = ct_fractions.loc[lambda x: x["cell_type"] == "Neutrophils", :]

# %%
datasets_with_neutros = (
    neutro_fractions.groupby("dataset")
    .apply(lambda x: np.sum(x["fraction"]) > 0)
    .where(lambda x: x)
    .dropna()
    .index.tolist()
)

# %%
neutro_subset = neutro_fractions.loc[
    lambda x: (
        x["dataset"].isin(datasets_with_neutros) & x["condition"].isin(["LUAD", "LSCC"])
    ),
    :,
].sort_values("fraction", ascending=False)
neutro_subset

# %%
fig, ax = plt.subplots()

neutro_subset["condition"] = neutro_subset["condition"].astype(str)
neutro_subset["dataset"] = neutro_subset["dataset"].astype(str)

sns.stripplot(data=neutro_subset, x="condition", y="fraction", hue="dataset", ax=ax, size=7, linewidth=2)
sns.boxplot(
    data=neutro_subset,
    x="condition",
    y="fraction",
    ax=ax,
    width=0.5
)
ax.legend(bbox_to_anchor=(1.1, 1.05))

# %%
mod = smf.ols("fraction ~ C(condition) + dataset", data=neutro_subset)
res = mod.fit()

# %%
res.summary()

# %%
res.pvalues["C(condition)[T.LUAD]"]

# %% [markdown]
# ## neutro recruitment signature
# (LSCC vs. LUAD) 

# %%
# TODO maybe best to also run DESeq2 for this

# %% [markdown]
# ---

# %%
adata.obs.groupby(
    ["dataset", "patient", "sample", "origin"], observed=True
).size().reset_index(name="n_cells").sort_values("n_cells", ascending=False).to_csv(
    "./granulocyte_count_per_patient.tsv", sep="\t"
)

# %%
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)

# %%
sc.tl.pca(adata, use_highly_variable=True)

# %%
sc.pp.neighbors(adata)

# %%
sc.tl.umap(adata)

# %%
sc.pl.umap(adata, color=["cell_type", "dataset"])

# %%
adata

# %%

# %% [markdown] tags=[]
# ---
#
# # Neutrophils

# %%
adata_neutro = adata[adata.obs["cell_type"] == "Granulocytes", :].copy()

# %%
sc.pp.neighbors(adata_neutro, use_rep="X_scANVI")
sc.tl.leiden(adata_neutro, resolution=0.5)
sc.tl.umap(adata_neutro)

# %%
sc.pl.umap(adata_neutro, color=["origin", "leiden", "dataset"], wspace=0.5)

# %%
sc.pl.umap(
    adata_neutro,
    color=adata.obs.columns[adata.obs.columns.str.startswith("scissor")],
    wspace=0.5,
    ncols=3,
)

# %%
neutro_ukimv = adata_neutro[adata_neutro.obs["dataset"] == "UKIM-V", :]

# %%
sc.pl.umap(
    neutro_ukimv,
    color=["origin", "leiden", "patient", "dataset", "VEGFA"],
    wspace=0.5,
    ncols=2,
)

# %%
sc.pl.dotplot(
    neutro_ukimv,
    groupby=["patient", "origin"],
    var_names="VEGFA",
    title="tumor_primary",
    vmin=0,
    vmax=1,
)

# %%
sc.pl.dotplot(
    neutro_ukimv[neutro_ukimv.obs["origin"] == "normal_adjacent", :],
    groupby="patient",
    var_names="VEGFA",
    title="normal_adjacent",
    vmin=0,
    vmax=1,
)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    sc.pl.umap(adata_neutro, color=["VEGFA"], cmap="inferno")

# %%
sc.pl.dotplot(
    adata,
    groupby=["patient", "origin"],
    var_names="VEGFA",
    title="tumor_primary",
    vmin=0,
    vmax=1,
)

# %%
adata_neutro.obs.groupby(
    ["dataset", "patient", "origin"], observed=True
).size().reset_index(name="n")

# %%
sc.pl.umap(adata_neutro, color=["origin", "dataset"])

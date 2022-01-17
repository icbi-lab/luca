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
adata_t = adata[adata.obs["cell_type"] == "Tumor cells", :].copy()
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
    "/home/sturm/Downloads/differential_signature_tumor_normal_progeny.tsv",
    sep="\t",
    index_col=0,
)
tumor_normal_cytosig = pd.read_csv(
    "/home/sturm/Downloads/differential_signature_tumor_normal_cytosig.tsv",
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
#
# All pvalues are unadjusted

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
tumor_vs_normal = ["PTGS2", "SELL", "CXCR2", "VEGFA", "OLR1", "CXCR4", "CXCR1", "ICAM1", "FCGR3B"]

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
deseq2_res_tumor_normal = pd.read_csv(
    "../../data/30_downstream_analyses/de_analysis/tumor_normal/de_deseq2/adata_tumor_normal_neutrophils_DESeq2_result.tsv",
    sep="\t",
    index_col=0,
).set_index("gene_id.1")

# %%
sc.pl.matrixplot(
    adata_tumor_normal, var_names=top_up_normal, groupby="origin", title="top up normal"
)

# %%
sc.pl.matrixplot(
    adata_tumor_normal, var_names=top_up_tumor, groupby="origin", title="top up tumor"
)

# %%
deseq2_res_tumor_normal

# %%
adata_primary = sh.pseudobulk.pseudobulk(
    adata[adata.obs["origin"] == "tumor_primary", :].copy(), groupby=["dataset", "patient", "cell_type_major"]
)
sc.pp.normalize_total(adata_primary, target_sum=1000)
sc.pp.log1p(adata_primary)

# %%
pd.set_option("display.max_rows", 200)

# %%
adata.var[adata.var_names.str.startswith("MT")]

# %%
adata.obs.columns

# %%
sc.pl.matrixplot(adata_primary, groupby="cell_type_major", var_names="VEGFA", swap_axes=True)

# %%
sh.pairwise.plot_paired(
    adata_tumor_normal,
    groupby="origin",
    paired_by="patient",
    var_names=tumor_vs_normal,
    pvalues=deseq2_res_tumor_normal.loc[tumor_vs_normal, "pvalue"],
    pvalue_template="DESeq2 p={:.2f}", 
    ylabel="log norm counts"
)

# %% [markdown]
# ### Using all samples (no pairing between tumor/normal) 
#
# This is a use-case of a linear mixed effects model

# %%
# tmp_adata = adata_n[
#     adata_n.obs["origin"].isin(["normal", "tumor_primary"]),
#     :,
# ]
# tmp_adata.obs["origin"] = [
#     "normal" if "normal" in x else x for x in tmp_adata.obs["origin"]
# ]
# adata_tumor_normal = sh.pseudobulk.pseudobulk(
#     tmp_adata, groupby=["dataset", "patient", "origin"]
# )
# sc.pp.normalize_total(adata_tumor_normal, target_sum=1000)
# sc.pp.log1p(adata_tumor_normal)

# %%
me_data = adata_tumor_normal.obs.join(pd.DataFrame(adata_tumor_normal.X, columns=adata_tumor_normal.var_names, index=adata_tumor_normal.obs_names))

# %%
me_pvalues = []
for gene in tumor_vs_normal:
    mod = smf.mixedlm(f"{gene} ~ origin", me_data, groups=me_data["dataset"])
    res = mod.fit()
    me_pvalues.append(res.pvalues["origin[T.tumor_primary]"])

# %%
sh.pairwise.plot_paired(
    adata_tumor_normal,
    groupby="origin",
    var_names=tumor_vs_normal,
    hue="dataset",
    show_legend=False,
    size=5,
    pvalues=me_pvalues,
    ylabel="log norm counts",
    pvalue_template="LME p={:.3f}"
)

# %% [markdown]
# ## Top DE genes as determined by DESeq2

# %%
sc.pl.matrixplot(
    adata_tumor_normal,
    var_names=deseq2_res_tumor_normal.index.values[:20],
    groupby="origin",
    title="top DE genes",
)

# %%
tmp_top10 = deseq2_res_tumor_normal.index.values[:10]
sh.pairwise.plot_paired(
    adata_tumor_normal,
    groupby="origin",
    paired_by="patient",
    var_names=tmp_top10,
    pvalues=deseq2_res_tumor_normal.loc[tmp_top10, "pvalue"],
    pvalue_template="DESeq2 p={:.4f}",
    ylabel="log norm counts"
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

sns.stripplot(
    data=neutro_subset,
    x="condition",
    y="fraction",
    hue="dataset",
    ax=ax,
    size=7,
    linewidth=2,
)
sns.boxplot(
    data=neutro_subset, x="condition", y="fraction", ax=ax, width=0.5, fliersize=0
)
ax.legend(bbox_to_anchor=(1.1, 1.05))

# %% [markdown]
# ## T cell fractions by tumor type

# %%
tcell_fractions = ct_fractions.loc[lambda x: x["cell_type"].isin(["T cell CD8"]), :]

# %%
datasets_with_tcells = (
    tcell_fractions.groupby("dataset")
    .apply(lambda x: np.sum(x["fraction"]) > 0)
    .where(lambda x: x)
    .dropna()
    .index.tolist()
)

# %%
tcell_subset = tcell_fractions.loc[
    lambda x: (
        x["dataset"].isin(datasets_with_tcells) & x["condition"].isin(["LUAD", "LSCC"])
    ),
    :,
].sort_values("fraction", ascending=False)
tcell_subset.head()

# %%
fig, ax = plt.subplots()

tcell_subset["condition"] = tcell_subset["condition"].astype(str)
tcell_subset["dataset"] = tcell_subset["dataset"].astype(str)

sns.stripplot(
    data=tcell_subset,
    x="condition",
    y="fraction",
    hue="dataset",
    ax=ax,
    size=7,
    linewidth=2,
)
sns.boxplot(
    data=tcell_subset, x="condition", y="fraction", ax=ax, width=0.5, fliersize=0
)
ax.legend(bbox_to_anchor=(1.1, 1.05))

# %%
mod = smf.ols("fraction ~ C(condition) + dataset", data=tcell_subset)
res = mod.fit()

# %%
res.summary()

# %%
res.pvalues["C(condition)[T.LUAD]"]

# %% [markdown]
# ## neutro recruitment signature
# (LSCC vs. LUAD) 

# %%
recruitment_genes = [
    "SOX2",
    "NKX2-1",
    "CXCL1",
    "CXCL2",
    "CXCL5",
    "CXCL6",
    "CXCL8",
    "IL6",
    # "IL17A", (<10 reads)
    "TNF",
    "LTA",
    "CSF2",
    "CSF3",
    "CCL2",
    "CCL3",
    "CCL4",
    "CCL5",
    "CXCL12",
]

# %%
deseq2_res_luad_lscc = pd.read_csv(
    "../../data/30_downstream_analyses/de_analysis/luad_lscc/de_deseq2/adata_luad_lscc_tumor_cells_DESeq2_result.tsv",
    sep="\t",
    index_col=0,
).set_index("gene_id.1")

# %%
tumor_cells_by_origin = sh.pseudobulk.pseudobulk(
    adata_t[
        (adata_t.obs["origin"] == "tumor_primary")
        & adata_t.obs["condition"].isin(["LUAD", "LSCC"]),
        :,
    ].copy(),
    groupby=["patient", "condition", "dataset"],
)
sc.pp.normalize_total(tumor_cells_by_origin, target_sum=1000)
sc.pp.log1p(tumor_cells_by_origin)

# %%
sc.pl.matrixplot(
    tumor_cells_by_origin, var_names=recruitment_genes, groupby="condition"
)

# %%
sh.pairwise.plot_paired(
    tumor_cells_by_origin,
    groupby="condition",
    var_names=recruitment_genes,
    hue="dataset",
    show_legend=False,
    size=5,
    ylabel="log norm counts",
    pvalues=deseq2_res_luad_lscc.loc[recruitment_genes, "padj"],
    pvalue_template="DESeq2 FDR={:.3f}"
)

# %%

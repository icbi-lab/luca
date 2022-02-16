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
import matplotlib
import altair as alt

# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
adata = sc.read_h5ad(
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad"
)

# %%
adata.obs.columns

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true
# # mRNA content
#
# Need to compute ratios, because the baseline difference between datasets and platforms is very high. 

# %%
rel_counts = (
    adata.obs.groupby(["dataset", "cell_type_coarse"])
    .agg(total_counts=("total_counts", "median"))
    .reset_index()
    .groupby("dataset")
    .apply(
        lambda x: x.assign(
            rel_counts=np.log2(x["total_counts"])
            - np.log2(
                x.loc[x["cell_type_coarse"] == "Epithelial cell", "total_counts"].values
            )
        )
    )
)

# %%
rel_counts

# %%
order = (
    rel_counts.groupby("cell_type_coarse")
    .mean()
    .sort_values("rel_counts")
    .index.tolist()
)
(
    alt.Chart(
        rel_counts,
        title="Mean detected counts per cell-type, relative to Epithelial cells",
    )
    .mark_bar()
    .encode()
    .encode(
        x=alt.X("cell_type_coarse", sort=order),
        y=alt.Y("mean(rel_counts):Q", title="log2(ratio)"),
        color=alt.Color(
            "cell_type_coarse",
            scale=sh.colors.altair_scale("cell_type_coarse"),
            legend=None,
        ),
    )
    + alt.Chart(rel_counts)
    .mark_errorbar(extent="ci")
    .encode(
        x=alt.X("cell_type_coarse", sort=order),
        y=alt.Y("rel_counts", title="log2(ratio)"),
    )
)

# %%
adata.obs.columns

# %%
for title, tmp_adata in {
    "counts per platform": adata,
    "counts per platform (Epithelial cells)": adata[
        adata.obs["cell_type_coarse"] == "Epithelial cell", :
    ],
}.items():
    counts_per_platform = (
        tmp_adata.obs.groupby(["sample", "platform"], observed=True)
        .agg(total_counts=("total_counts", "median"))
        .reset_index()
    )
    order = (
        counts_per_platform.groupby("platform")
        .median()
        .sort_values("total_counts")
        .index.tolist()
    )
    alt.Chart(counts_per_platform, title=title).mark_boxplot().encode(
        x=alt.X("platform", sort=order[::-1]),
        y=alt.Y("total_counts", scale=alt.Scale(type="log")),
        color=alt.Color(
            "platform",
            scale=sh.colors.altair_scale("platform"),
            legend=None,
        ),
    ).display()

# %% [markdown]
# # Neutrophil subset

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
# ## UKIM-V datasets only

# %%
adata_n_ukimv = adata_n[adata_n.obs["dataset"].str.startswith("UKIM-V"), :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_n_ukimv, use_rep="X_scANVI")

# %%
adata_n_ukimv.obs

# %%
with plt.rc_context({"figure.figsize": (4, 4), "figure.dpi": 300}):
    sc.pl.umap(adata_n_ukimv, color="origin", frameon=False, size=20)

# %% [markdown]
# # DE analysis tumor normal

# %%
# candidate genes
tumor_vs_normal = [
    "PTGS2",
    "SELL",
    "CXCR2",
    "VEGFA",
    "OLR1",
    "CXCR4",
    "CXCR1",
    "ICAM1",
    "FCGR3B",
    "CD83",
    "ARG1",
    "CCL2",
    "JUN",
]

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
sc.pp.normalize_total(adata_tumor_normal, target_sum=1e6)
adata_tumor_normal.layers["cpm"] = adata_tumor_normal.X.copy()
sc.pp.log1p(adata_tumor_normal)
adata_tumor_normal.obs["origin2"] = [
    "T" if o == "tumor_primary" else "N" for o in adata_tumor_normal.obs["origin"]
]

# %% [markdown]
# ### Overview pseudobulk samples

# %%
adata_tumor_normal.obs.loc[
    lambda x: x["dataset"].str.startswith("UKIM"), :
].sort_values(["dataset", "patient", "origin"])

# %%
deseq2_res_tumor_normal = pd.read_csv(
    "../../data/30_downstream_analyses/de_analysis/tumor_normal/de_deseq2/adata_tumor_normal_neutrophils_DESeq2_result.tsv",
    sep="\t",
    index_col=0,
).set_index("gene_id.1")

# %%
deseq2_res_tumor_normal

# %%
sh.pairwise.plot_paired(
    adata_tumor_normal,
    groupby="origin2",
    paired_by="patient",
    var_names=tumor_vs_normal,
    pvalues=deseq2_res_tumor_normal.loc[tumor_vs_normal, "padj"],
    pvalue_template="DESeq2 FDR={:.2f}",
    ylabel="log norm counts",
    n_cols=5,
)

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# ### Using all samples (no pairing between tumor/normal) 
#
# This is a use-case of a linear mixed effects model

# %%
me_data = adata_tumor_normal.obs.join(
    pd.DataFrame(
        adata_tumor_normal.X,
        columns=adata_tumor_normal.var_names,
        index=adata_tumor_normal.obs_names,
    )
)

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
    pvalue_template="LME p={:.3f}",
)

# %% [markdown]
# ## Top DE genes as determined by DESeq2

# %%
tmp_top10 = deseq2_res_tumor_normal.sort_values("padj").index.values[:10]
sh.pairwise.plot_paired(
    adata_tumor_normal,
    groupby="origin2",
    paired_by="patient",
    var_names=tmp_top10,
    pvalues=deseq2_res_tumor_normal.loc[tmp_top10, "padj"],
    pvalue_template="DESeq2 FDR={:.4f}",
    ylabel="log norm counts",
)

# %%
tmp_top = deseq2_res_tumor_normal.sort_values("padj").index.values[:30]
sh.pairwise.plot_paired_fc(
    adata_tumor_normal,
    groupby="origin2",
    paired_by="patient",
    var_names=tmp_top,
    layer="cpm",
).properties(height=150)

# %% [markdown]
# # neutro fractions by tumor type

# %%
ct_fractions = (
    adata[
        (adata.obs["origin"] == "tumor_primary")
        & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"]),
        :,
    ]
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
ax.set_title("Neutrophil fraction in LSCC vs LUAD")

# %%
mod = smf.ols("fraction ~ C(condition) + dataset", data=neutro_subset)
res = mod.fit()

# %%
res.pvalues

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# # T cell fractions by tumor type

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
# ## neutro recruitment signature (LSCC vs. LUAD) 

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

# %% [markdown]
# ### genes of interest

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
    pvalue_template="DESeq2 FDR={:.3f}",
    n_cols=5,
)

# %% [markdown]
# ### top genes

# %%
top_genes = deseq2_res_luad_lscc.index[:30]
sh.pairwise.plot_paired(
    tumor_cells_by_origin,
    groupby="condition",
    var_names=top_genes,
    hue="dataset",
    show_legend=False,
    size=5,
    ylabel="log norm counts",
    pvalues=deseq2_res_luad_lscc.loc[top_genes, "padj"],
    pvalue_template="DESeq2 FDR={:.3f}",
    n_cols=10,
)

# %%
sox2_genes = [
    "ABCC6",
    "CCND1",
    "DLGAP1",
    "GLI2",
    "GLI3",
    "HHAT",
    # 'ISL1',
    "NANOG",
    "NTRK3",
    # "RHO",
]
sh.pairwise.plot_paired(
    tumor_cells_by_origin,
    groupby="condition",
    var_names=sox2_genes,
    hue="dataset",
    show_legend=False,
    size=5,
    ylabel="log norm counts",
    pvalues=deseq2_res_luad_lscc.loc[sox2_genes, "padj"],
    pvalue_template="DESeq2 FDR={:.3f}",
    n_cols=5,
)

# %% [markdown]
# # VEGFA sources in NSCLC

# %%
adata_primary = sh.pseudobulk.pseudobulk(
    adata[
        (adata.obs["origin"] == "tumor_primary")
        & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"]),
        :,
    ].copy(),
    groupby=["dataset", "patient", "cell_type_major"],
    min_obs=20,
)
sc.pp.normalize_total(adata_primary, target_sum=1e6)
sc.pp.log1p(adata_primary)

# %%
pd.set_option("display.max_rows", 200)

# %%
adata.obs.columns

# %%
df = adata_primary.obs
df["VEGFA"] = adata_primary[:, "VEGFA"].X

# %%
df

# %%
order = (
    df.groupby("cell_type_major")
    .agg("median")
    .sort_values("VEGFA", ascending=False)
    .index.values
)

# %%
PROPS = {
    "boxprops": {"facecolor": "none", "edgecolor": "black"},
    "medianprops": {"color": "black"},
    "whiskerprops": {"color": "black"},
    "capprops": {"color": "black"},
}

fig, ax = plt.subplots(1, 1, figsize=(10, 5))
sns.stripplot(
    x="cell_type_major", y="VEGFA", hue="dataset", data=df, ax=ax, order=order
)
sns.boxplot(
    x="cell_type_major",
    y="VEGFA",
    ax=ax,
    data=df,
    order=order,
    color="white",
    **PROPS,
    showfliers=False
)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
_ = plt.xticks(rotation=90)

# %% [markdown]
# # Compare cytosig and progeny

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

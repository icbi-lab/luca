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
from pathlib import Path
from tqdm.auto import tqdm

# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %%
signatures = pd.read_csv(
    "../../tables/gene_annotations/neutro_signatures.csv", comment="#"
).loc[lambda x: x["signature"] != "sig_pd1"]

# %%
ah = AnnotationHelper()

# %%
adata_n = sc.read_h5ad(
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/adata_neutrophil_clusters.h5ad"
)

# %%
immunomodulatory_genes = pd.read_excel("../../tables/gene_annotations/immunomodulatory_genes.xlsx")

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
autophagy_genes = """ELANE
CTSG
ITGB2
LCN2
ADAMTS12
ADAMTS13
ADAM10
ADAM11
ADAM12
ADAM17
ADAMTS4
ADAM28
LOXL2
MMP2
MMP3
ILF3
MMP7
MMP8
MMP9
MMP10
MMP11
MMP12
MMP13
MMP14
MMP15
MMP16
MMP17
TIMP1
TIMP2
TIMP3
TIMP4
""".split()

# %% [markdown]
# # UMAPs by covariate

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(adata_n, color="cell_type", legend_loc="on data", legend_fontoutline=2, frameon=False, size=20)

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(adata_n, color="condition", legend_fontoutline=2, frameon=False, size=20)

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(adata_n[adata_n.obs["origin"] != "nan", :], color="origin", legend_fontoutline=2, frameon=False, size=20)

# %%
adata_n.obs["dataset"].astype(str).unique()

# %%
adata_n.obs["origin_biopsy"] = adata_n.obs["origin"].astype(str)
adata_n.obs.loc[adata_n.obs["dataset"].isin(["Wu_Zhou_2021", "Kim_Lee_2020"]) & (adata_n.obs["origin"] == "tumor_primary"), "origin_biopsy"] = "biopsy"

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(adata_n, color="origin_biopsy", legend_fontoutline=2, frameon=False, size=20, groups="biopsy")

# %% [markdown]
# # Clusters by patient and dataset

# %%
patient_fracs = (
    adata_n.obs.groupby(["cell_type"])
    .apply(lambda x: x["patient"].value_counts(normalize=True))
    .unstack()
    .melt(ignore_index=False, var_name="patient", value_name="fraction")
    .reset_index()
)

dataset_fracs = (
    adata_n.obs.groupby(["cell_type"])
    .apply(lambda x: x["dataset"].value_counts(normalize=True))
    .unstack()
    .melt(ignore_index=False, var_name="dataset", value_name="fraction")
    .reset_index()
)

# %%
alt.Chart(dataset_fracs).mark_bar().encode(x="fraction", y="cell_type", color="dataset")

# %%
alt.Chart(patient_fracs).mark_bar().encode(x=alt.X("fraction", scale=alt.Scale(domain=[0,1])), y="cell_type", color="patient")

# %% [markdown]
# # UMAP by candidate genes

# %%
pb_n = sh.pseudobulk.pseudobulk(adata_n, groupby=["cell_type", "patient"])

# %%
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n)

# %%
sc.pl.umap(
    adata_n,
    color=tumor_vs_normal,
    cmap="inferno",
    size=20,
    ncols=5,
    frameon=False,
)

# %% [markdown]
# ### genes of interest (2)

# %%
sc.pl.umap(
    adata_n,
    color=["OLR1", "CD36", "ITGA2B", "ITGB3"],
    cmap="inferno",
    size=20,
    ncols=5,
    frameon=False,
)

# %%
sc.pl.dotplot(
    adata_n,
    var_names=["OLR1", "CD36", "ITGA2B", "ITGB3"],
    groupby=["cell_type"],
)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=["OLR1", "CD36", "ITGA2B", "ITGB3"],
    groupby=["cell_type"],
    cmap="bwr"
)

# %% [markdown]
# ### genes of interest (3) -- autophagy

# %%
sc.pl.dotplot(
    adata_n,
    var_names=autophagy_genes,
    groupby=["cell_type"],
)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=autophagy_genes,
    groupby=["cell_type"],
    cmap="bwr"
)

# %% [markdown]
# ## Genes of interest (4) - Immunomodulatory

# %%
sc.pl.dotplot(
    adata_n,
    var_names=immunomodulatory_genes.loc[lambda x: x["type"] == "immunomodulatory", "gene_symbol"],
    groupby=["cell_type"]
)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=immunomodulatory_genes.loc[lambda x: x["type"] == "immunomodulatory", "gene_symbol"],
    groupby=["cell_type"], 
    cmap="bwr"
)

# %% [markdown]
# ## Genes of interest (5) - Immune response

# %%
sc.pl.dotplot(
    adata_n,
    var_names=immunomodulatory_genes.loc[lambda x: x["type"] == "immune_response", "gene_symbol"],
    groupby=["cell_type"]
)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=immunomodulatory_genes.loc[lambda x: x["type"] == "immune_response", "gene_symbol"],
    groupby=["cell_type"],
    cmap="bwr"
)

# %% [markdown]
# ---

# %% [markdown]
# # Find marker genes for Neutrophil clusters

# %%
sc.tl.rank_genes_groups(pb_n, groupby="cell_type", method="t-test", use_raw=False)

# %%
sc.tl.rank_genes_groups(adata_n, "cell_type")

# %%
sc.pl.rank_genes_groups_matrixplot(pb_n, dendrogram=False, cmap="bwr")

# %%
sc.pl.rank_genes_groups_dotplot(adata_n, dendrogram=False)

# %%
sc.pl.umap(
    adata_n,
    color=["CXCR2", "S100A12", "CTSC", "CXCL2", "IFIT1", "RPL5"],
    cmap="inferno",
    size=20,
    ncols=3,
    frameon=False,
)

# %%
sc.pl.dotplot(
    adata_n,
    var_names=["CXCR2", "S100A12", "CTSC", "CXCL2", "IFIT1", "RPL5"],
    groupby="cell_type",
)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=["CXCR2", "S100A12", "CTSC", "CXCL2", "IFIT1", "RPL5"],
    groupby="cell_type",
    cmap="bwr",
)

# %%
signature_dict = {}
for sig in signatures["signature"].unique():
    signature_dict[sig] = signatures.loc[
        lambda x: x["signature"] == sig, "gene_symbol"
    ].tolist()

# %%
sc.pl.matrixplot(pb_n, var_names=signature_dict, cmap="bwr", groupby="cell_type")

# %%
sc.pl.matrixplot(
    pb_n, var_names=signature_dict["tan_sig"], cmap="bwr", groupby="cell_type"
)

# %%
for sig, genes in signature_dict.items():
    sc.tl.score_genes(pb_n, genes, score_name=sig)
    sc.tl.score_genes(adata_n, genes, score_name=sig)

# %%
sc.pl.matrixplot(
    pb_n, var_names=list(signature_dict.keys()), cmap="bwr", groupby="cell_type"
)

# %%
sc.pl.matrixplot(
    adata_n[adata_n.obs["cell_type"] != "other"],
    var_names=list(signature_dict.keys()),
    cmap="bwr",
    groupby="cell_type",
)

# %%
sc.pl.umap(adata_n, color=list(signature_dict.keys()), cmap="inferno", size=20, ncols=3)

# %% [markdown]
# ## UKIM-V datasets only

# %%
adata_n_ukimv = adata_n[adata_n.obs["dataset"].str.startswith("UKIM-V"), :].copy()

# %%
ah.reprocess_adata_subset_scvi(adata_n_ukimv, use_rep="X_scANVI")

# %%
adata_n_ukimv.obs

# %%
tmp_df = (
    adata_n.obs.groupby(["cell_type"])
    .apply(lambda x: x["origin"].value_counts(normalize=True))
    .unstack()
    .melt(var_name="origin", value_name="fraction", ignore_index=False)
    .reset_index()
    .loc[lambda x: x["cell_type"] != "other"]
)

# %%
alt.Chart(tmp_df).mark_bar().encode(x="fraction", y="cell_type", color="origin")

# %%
tmp_df = (
    adata_n.obs.groupby(["cell_type"])
    .apply(lambda x: x["condition"].value_counts(normalize=True))
    .unstack()
    .melt(var_name="condition", value_name="fraction", ignore_index=False)
    .reset_index()
    .loc[lambda x: x["cell_type"] != "other"]
)

# %%
alt.Chart(tmp_df).mark_bar().encode(x="fraction", y="cell_type", color="condition")

# %% [markdown]
# # DE analysis tumor normal

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

# %%
sc.pl.umap(
    adata_n,
    color=deseq2_res_tumor_normal.sort_values("padj").index.values[:30],
    cmap="inferno",
    size=20,
    ncols=5,
)

# %% [markdown]
# # LUAD vs. LUSC
#
# ## Neutro fractions

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
        x["dataset"].isin(datasets_with_neutros) & x["condition"].isin(["LUAD", "LUSC"])
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
ax.set_title("Neutrophil fraction in LUSC vs LUAD")

# %%
mod = smf.ols("fraction ~ C(condition) + dataset", data=neutro_subset)
res = mod.fit()

# %%
res.pvalues

# %% [markdown]
# ## neutro recruitment signature (LUSC vs. LUAD)

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
    "CSF1",
    "CSF2",
    "CSF3",
    "CCL2",
    "CCL3",
    "CCL4",
    "CCL5",
    "CXCL12",
]

# %%
deseq2_res_luad_lusc = pd.read_csv(
    "../../data/30_downstream_analyses/de_analysis/luad_lusc/de_deseq2/adata_primary_tumor_tumor_cells_DESeq2_result.tsv",
    sep="\t",
    index_col=0,
).set_index("gene_id.1")

# %%
tumor_cells_by_origin = sh.pseudobulk.pseudobulk(
    adata_t[
        (adata_t.obs["origin"] == "tumor_primary")
        & adata_t.obs["condition"].isin(["LUAD", "LUSC"]),
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
    pvalues=deseq2_res_luad_lusc.loc[recruitment_genes, "padj"],
    pvalue_template="DESeq2 FDR={:.3f}",
    n_cols=5,
)

# %% [markdown]
# ### top genes

# %%
top_genes = deseq2_res_luad_lusc.index[:30]
sh.pairwise.plot_paired(
    tumor_cells_by_origin,
    groupby="condition",
    var_names=top_genes,
    hue="dataset",
    show_legend=False,
    size=5,
    ylabel="log norm counts",
    pvalues=deseq2_res_luad_lusc.loc[top_genes, "padj"],
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
    pvalues=deseq2_res_luad_lusc.loc[sox2_genes, "padj"],
    pvalue_template="DESeq2 FDR={:.3f}",
    n_cols=5,
)

# %% [markdown]
# ### Signalling

# %%
interactions_of_interest = pd.read_excel(
    "../../tables/gene_annotations/neutro_recruitment_chemokines.xlsx"
)

# %%
dfs = []
for p in Path(
    "../../data/30_downstream_analyses/de_analysis/luad_lusc/de_deseq2/"
).glob("*.tsv"):
    cell_type = p.stem.replace("adata_primary_tumor_", "").replace("_DESeq2_result", "")
    dfs.append(pd.read_csv(p, sep="\t").assign(cell_type=cell_type))

# %%
selected_de_res = (
    pd.concat(dfs)
    .drop("gene_id", axis="columns")
    .rename(columns={"gene_id.1": "gene"})
    .loc[lambda x: x["gene"].isin(interactions_of_interest["tme_ligand"])]
    .loc[lambda x: ~pd.isnull(x["pvalue"])]
    .pipe(sh.util.fdr_correction)
    .sort_values("fdr")
    .copy()
)

# %%
selected_de_res.pipe(
    sh.compare_groups.pl.plot_lm_result_altair,
    x="gene",
    y="cell_type",
    color="log2FoldChange",
)

# %% [markdown]
# # early vs late

# %%
adata_primary = sh.pseudobulk.pseudobulk(
    adata[
        (adata.obs["origin"] == "tumor_primary")
        & ~adata.obs["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"]),
        :,
    ].copy(),
    groupby=["dataset", "patient", "condition", "tumor_stage", "cell_type_major"],
    min_obs=20,
)
sc.pp.normalize_total(adata_primary, target_sum=1e6)
sc.pp.log1p(adata_primary)

# %%
sh.pairwise.plot_paired(
    adata_primary[adata_primary.obs["cell_type_major"] == "Neutrophils", :],
    groupby="tumor_stage",
    var_names=["IGHA2", "DDX3Y"],
    ylabel="log norm counts",
    hue="dataset",
    panel_size=(7, 4),
)

# %% [markdown]
# # VEGFA sources in NSCLC

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
    showfliers=False,
)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
_ = plt.xticks(rotation=90)

# %% [markdown]
# # Compare cytosig and progeny

# %%
tumor_normal_progeny = pd.read_csv(
    "../../data/30_downstream_analyses/plots_and_comparisons/91_compare_groups/artifacts/tumor_normal_progeny.tsv",
    sep="\t",
    index_col=0,
)
tumor_normal_cytosig = pd.read_csv(
    "../../data/30_downstream_analyses/plots_and_comparisons/91_compare_groups/artifacts/tumor_normal_cytosig.tsv",
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
# # Differences in TANs between LUAD and LUSC

# %%
neutrophils_by_origin = sh.pseudobulk.pseudobulk(
    adata_n[
        (adata_n.obs["origin"] == "tumor_primary")
        & adata_n.obs["condition"].isin(["LUAD", "LUSC"]),
        :,
    ].copy(),
    groupby=["patient", "condition", "dataset"],
)
sc.pp.normalize_total(neutrophils_by_origin, target_sum=1000)
sc.pp.log1p(neutrophils_by_origin)

# %%
neutrophils_by_origin.obs.sort_values(["dataset", "condition"])

# %%
deseq2_neutro_luad_lusc = pd.read_csv(
    "../../data/30_downstream_analyses/de_analysis/luad_lusc/de_deseq2/adata_luad_lusc_neutrophils_DESeq2_result.tsv",
    sep="\t",
    index_col=0,
).sort_values("padj")

# %%
top_genes = deseq2_neutro_luad_lusc["gene_id.1"][:30].values

# %%
top_genes = deseq2_neutro_luad_lusc["gene_id.1"][:30]
sh.pairwise.plot_paired(
    neutrophils_by_origin,
    groupby="condition",
    var_names=top_genes,
    hue="dataset",
    show_legend=False,
    size=5,
    ylabel="log norm counts",
    pvalues=deseq2_neutro_luad_lusc.set_index("gene_id.1").loc[top_genes, "padj"],
    pvalue_template="DESeq2 FDR={:.3f}",
    n_cols=10,
)

# %%
deseq2_neutro_luad_lusc.iloc[:30].rename(columns={"gene_id.1": "gene"})

# %%
df = (
    deseq2_neutro_luad_lusc.loc[lambda x: x["padj"] < 0.1]
    .rename(columns={"gene_id.1": "gene"})
    .assign(
        ymin=lambda x: x["log2FoldChange"] - x["lfcSE"],
        ymax=lambda x: x["log2FoldChange"] + x["lfcSE"],
    )
)
order = df.sort_values("log2FoldChange")["gene"].values.tolist()[::-1]
(
    alt.Chart(df)
    .mark_bar()
    .encode(
        x=alt.X("gene", sort=order),
        y=alt.Y("log2FoldChange", title="log2 fold change"),
        color=alt.Color(
            "log2FoldChange", scale=alt.Scale(scheme="redblue", reverse=True)
        ),
    )
    .properties(height=100)
) + alt.Chart(df).mark_errorbar().encode(
    x=alt.X("gene", sort=order),
    y=alt.Y("ymin:Q", title="log2 fold change"),
    y2="ymax:Q",
)

# %%

# %% [markdown]
# ---

# %%
signature_genes = []
for cl in tqdm(pb_n.obs["leiden"].unique()):
    signature_genes.append(
        pd.DataFrame().assign(
            gene_symbol=pb_n.var_names,
            fold_change=sh.signatures.fold_change(
                pb_n, obs_col="leiden", positive_class=cl, inplace=False
            ),
            specific_fold_change=sh.signatures.specific_fold_change(
                pb_n, obs_col="leiden", positive_class=cl, inplace=False
            ),
            roc_auc=sh.signatures.roc_auc(
                pb_n, obs_col="leiden", positive_class=cl, inplace=False
            ),
            leiden=cl,
        )
    )

# %%
signature_genes_df = (
    pd.concat(signature_genes)
    .loc[
        lambda x: (x["fold_change"] > 0.5)
        & (x["specific_fold_change"] > 0.5)
        & (x["roc_auc"] > 0.7)
    ]
    .sort_values(["roc_auc", "specific_fold_change"], ascending=False)
    .groupby("leiden")
    .apply(lambda x: x.head(10))
)

# %%
gene_dict_for_plot = {}
for leiden in signature_genes_df["leiden"].unique():
    gene_dict_for_plot[leiden] = signature_genes_df.loc[
        lambda x: x["leiden"] == leiden
    ]["gene_symbol"].tolist()

# %%
signature_genes_df

# %% [markdown]
# ---

# %% [markdown]
# ### Neutrophils

# %%
results["luad_lusc"]["cytosig"].loc[lambda x: x["cell_type"] == "Neutrophils", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Cytosig (Neutrophils cells)")

# %%
results["luad_lusc"]["dorothea"].loc[lambda x: x["cell_type"] == "Neutrophils", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="TFs (tumor cells LUAD/LUSC)")

# %%
results["early_advanced"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="TFs (Neutrophils tumor/normal)"
)

# %%
results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].to_csv("/home/sturm/Downloads/neutrophils_tfs.tsv", sep="\t")

# %%
results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="TFs (Neutrophils tumor/normal)"
)

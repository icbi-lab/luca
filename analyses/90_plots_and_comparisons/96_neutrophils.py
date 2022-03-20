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
import warnings
from tqdm.auto import tqdm
import itertools

# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %%
artifact_dir = "../../data/30_downstream_analyses/neutrophils"
cpus = 16

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
adata = sc.read_h5ad(
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/full_atlas_neutrophil_clusters.h5ad"
)

# %%
immunomodulatory_genes = pd.read_excel(
    "../../tables/gene_annotations/immunomodulatory_genes.xlsx"
)

# %%
neutro_recruitment_genes = pd.read_excel(
    "../../tables/gene_annotations/neutro_recruitment_chemokines.xlsx"
)

# %%

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

# %%
adata_n.obs["cell_type_tan_nan"] = [x[:4] for x in adata_n.obs["cell_type"]]

# %% [markdown]
# # UMAPs by covariate

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(
        adata_n,
        color="cell_type_tan_nan",
        legend_loc="on data",
        legend_fontoutline=2,
        frameon=False,
        size=20,
    )

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(
        adata_n,
        color="cell_type",
        legend_loc="on data",
        legend_fontoutline=2,
        frameon=False,
        size=20,
    )

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(adata_n, color="condition", legend_fontoutline=2, frameon=False, size=20)

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(
        adata_n[adata_n.obs["origin"] != "nan", :],
        color="origin",
        legend_fontoutline=2,
        frameon=False,
        size=20,
    )

# %%
adata_n.obs["dataset"].astype(str).unique()

# %%
adata_n.obs["origin_biopsy"] = adata_n.obs["origin"].astype(str)
adata_n.obs.loc[
    adata_n.obs["dataset"].isin(["Wu_Zhou_2021", "Kim_Lee_2020"])
    & (adata_n.obs["origin"] == "tumor_primary"),
    "origin_biopsy",
] = "biopsy"

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(
        adata_n,
        color="origin_biopsy",
        legend_fontoutline=2,
        frameon=False,
        size=20,
        groups="biopsy",
    )

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
alt.Chart(patient_fracs).mark_bar().encode(
    x=alt.X("fraction", scale=alt.Scale(domain=[0, 1])), y="cell_type", color="patient"
)

# %% [markdown]
# # Find marker genes for Neutrophil clusters

# %%
pb_n = sh.pseudobulk.pseudobulk(adata_n, groupby=["cell_type", "patient"])

# %%
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n, base=2)

# %%
pb_tan_nan = sh.pseudobulk.pseudobulk(adata_n, groupby=["cell_type_tan_nan", "patient"])

# %%
sc.pp.normalize_total(pb_tan_nan, target_sum=1e6)
sc.pp.log1p(pb_tan_nan, base=2)

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    marker_res = sh.compare_groups.lm.test_lm(
        pb_n,
        "~ C(cell_type, Sum) + patient",
        "cell_type",
    )

# %%
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    marker_res_tan_nan = sh.compare_groups.lm.test_lm(
        pb_tan_nan[pb_tan_nan.obs["patient"].duplicated(keep=False), :].copy(),
        "~ C(cell_type_tan_nan, Treatment('NAN-')) + patient",
        "cell_type_tan_nan",
        contrasts="Treatment('NAN-')",
    )

# %%
marker_res_tan_nan = (
    # the valuese are already on log-scale, therefore using logfun=identity
    marker_res_tan_nan.pipe(sh.util.log2_fc, logfun=lambda x: x)
    .loc[lambda x: ~pd.isnull(x["pvalue"])]
    .pipe(sh.util.fdr_correction)
    .sort_values("pvalue")
    .assign(group=lambda x: ["NAN" if _ > 0 else "TAN" for _ in x["log2_fc"]])
)

# %%
marker_res = (
    # the valuese are already on log-scale, therefore using logfun=identity
    marker_res.pipe(sh.util.log2_fc, logfun=lambda x: x)
    .loc[lambda x: ~pd.isnull(x["pvalue"]) & (x["log2_fc"] > 0)]
    .pipe(sh.util.fdr_correction)
    .sort_values("pvalue")
)

# %%
marker_res_tan_nan.to_csv(f"{artifact_dir}/tan_vs_nan_markers.csv")

# %%
marker_res.to_csv(f"{artifact_dir}/neutro_cluster_markers.csv")

# %% [markdown]
# ## Top markers TAN vs NAN

# %%
top_markers_tan_nan = (
    marker_res_tan_nan.loc[lambda x: (abs(x["log2_fc"]) > 1) & (x["fdr"] < 0.01)]
    .sort_values("log2_fc", key=abs, ascending=False)
    .groupby("group")
    .apply(lambda x: x.head(15))
)

# %%
sc.pl.matrixplot(
    pb_n,
    groupby="cell_type",
    var_names=top_markers_tan_nan["variable"][::-1],
    cmap="bwr",
)

# %%
sh.pairwise.plot_paired_fc(
    pb_tan_nan,
    "cell_type_tan_nan",
    paired_by="patient",
    var_names=top_markers_tan_nan.iloc[:30].sort_values("log2_fc")["variable"],
    metric="diff",  # diff metric = log fold change, as data is already log-transformed.
    metric_name="log2 fold change",
)

# %% [markdown]
# ## Top markers Neutrophil subclusters

# %%
top_markers = (
    marker_res.loc[lambda x: (x["fdr"] < 0.01) & (x["log2_fc"] > 1) & (x["coef"] > 1)]
    .sort_values("log2_fc", ascending=False)
    .groupby("group")
    .apply(lambda x: x.head(15))
)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names={
        g: top_markers.loc[lambda x: x["group"] == g, "variable"]
        for g in top_markers["group"].unique()
    },
    groupby="cell_type",
    cmap="bwr",
)

# %% [markdown]
# ## Selected Top markers

# %%
selected_top_markers = {
    "NAN": ["CST7", "TMX4"],
    "TAN": ["CCL4", "CCL3"],
    "NAN-1": ["S100A12", "RBP7"],
    "NAN-2": ["ENTPD1", "SLC8A1"],
    "TAN-1": ["PIGR", "CD52"],
    "TAN-2": ["PRR5L", "IFIT1"],
    "TAN-3": ["PLPP3", "CD22"],
}

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=selected_top_markers,
    groupby="cell_type",
    cmap="bwr",
)

# %%
sc.pl.umap(
    adata_n,
    color=itertools.chain.from_iterable(selected_top_markers.values()),
    cmap="inferno",
    size=40,
    ncols=5,
    frameon=False,
)

# %% [markdown]
# ## Pathways, TFs, CytoSig

# %%
res = sh.compare_groups._run_tool(adata_n, sh.compare_groups.TOOLS)

# %%
results = sh.compare_groups.compare_signatures(
    "neutrophil_clusters",
    {tool: {"neutrophils": ad} for tool, ad in res.items()},
    n_jobs=cpus,
    pseudobulk_group_by=["patient"],
    column_to_test="cell_type",
    contrasts="Sum",
    lm_covariate_str="+ patient",
)

# %%
sh.compare_groups.pl.plot_lm_result_altair(results["progeny"], p_cutoff=1, title="Progeny Neutrophils")

# %%
sh.compare_groups.pl.plot_lm_result_altair(results["dorothea"], title="Dorothea Neutrophils", p_cutoff=0.01)

# %%
sh.compare_groups.pl.plot_lm_result_altair(results["cytosig"], title="Cytosig Neutrophils", p_cutoff=0.01)

# %% [markdown] tags=[]
# # UMAP and Heatmaps of selected gene sets

# %% [markdown]
# ## TAN vs NAN selected genes

# %%
sc.pl.umap(
    adata_n,
    color=tumor_vs_normal,
    cmap="inferno",
    size=20,
    ncols=5,
    frameon=False,
)

# %%
sc.pl.matrixplot(pb_n, var_names=tumor_vs_normal, groupby=["cell_type"], cmap="bwr")

# %% [markdown]
# ## genes of interest (2)

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
    cmap="bwr",
)

# %% [markdown]
# ## genes of interest (3) -- autophagy

# %%
sc.pl.dotplot(
    adata_n,
    var_names=autophagy_genes,
    groupby=["cell_type"],
)

# %%
sc.pl.matrixplot(pb_n, var_names=autophagy_genes, groupby=["cell_type"], cmap="bwr")

# %% [markdown]
# ## Genes of interest (4) - Immunomodulatory

# %%
sc.pl.dotplot(
    adata_n,
    var_names=immunomodulatory_genes.loc[
        lambda x: x["type"] == "immunomodulatory", "gene_symbol"
    ],
    groupby=["cell_type"],
)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=immunomodulatory_genes.loc[
        lambda x: x["type"] == "immunomodulatory", "gene_symbol"
    ],
    groupby=["cell_type"],
    cmap="bwr",
)

# %% [markdown]
# ## Genes of interest (5) - Immune response

# %%
sc.pl.dotplot(
    adata_n,
    var_names=immunomodulatory_genes.loc[
        lambda x: x["type"] == "immune_response", "gene_symbol"
    ],
    groupby=["cell_type"],
)

# %%
sc.pl.matrixplot(
    pb_n,
    var_names=immunomodulatory_genes.loc[
        lambda x: x["type"] == "immune_response", "gene_symbol"
    ],
    groupby=["cell_type"],
    cmap="bwr",
)

# %% [markdown]
# # LUAD vs. LUSC
#
# ## Neutro fractions

# %%
ct_fractions = (
    adata[
        (adata.obs["origin"] == "tumor_primary")
        & ~adata.obs["dataset"].isin(["Guo_Zhang_2018"]),
        :,
    ]
    .obs.groupby(["dataset", "patient", "condition", "study"])["cell_type_major"]
    .value_counts(normalize=True)
    .reset_index()
)

# %%
ct_fractions = ct_fractions.rename(
    columns={"level_4": "cell_type_major", "cell_type_major": "fraction"}
)

# %%
neutro_fractions = ct_fractions.loc[lambda x: x["cell_type_major"] == "Neutrophils", :]

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
mod = smf.ols("fraction ~ C(condition) + dataset", data=neutro_subset)
res = mod.fit()

# %%
fig, ax = plt.subplots()

neutro_subset["condition"] = neutro_subset["condition"].astype(str)
neutro_subset["dataset"] = neutro_subset["dataset"].astype(str)

sns.stripplot(
    data=neutro_subset,
    x="condition",
    y="fraction",
    hue="study",
    palette=sh.colors.COLORS.study,
    ax=ax,
    size=7,
    linewidth=2,
)
sns.boxplot(
    data=neutro_subset, x="condition", y="fraction", ax=ax, width=0.5, fliersize=0
)
ax.legend(bbox_to_anchor=(1.1, 1.05))
ax.set_title("Neutrophil fraction in LUSC vs LUAD\np={:.4f}, linear model".format(res.pvalues["C(condition)[T.LUSC]"]))
plt.show()

# %% [markdown]
# ## neutro recruitment signature (LUSC vs. LUAD)

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
# # VEGFA sources in NSCLC

# %%
pb_primary = sh.pseudobulk.pseudobulk(
    adata[adata.obs["origin"] == "tumor_primary", :].copy(),
    groupby=["cell_type_major", "patient", "study"],
)

# %%
sc.pp.normalize_total(pb_primary, target_sum=1e6)
sc.pp.log1p(pb_primary)

# %%
df = pb_primary.obs
df["VEGFA"] = pb_primary[:, "VEGFA"].X

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
    x="cell_type_major",
    y="VEGFA",
    hue="study",
    data=df,
    ax=ax,
    order=order,
    palette=sh.colors.COLORS.study,
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
# # Signature scores

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

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
from IPython.display import display_html
import itertools

# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %%
artifact_dir = "../../data/30_downstream_analyses/neutrophils"
cpus = 16

# %%
ah = AnnotationHelper()

# %%
adata_n = sc.read_h5ad(
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/adata_neutrophil_clusters.h5ad"
    # "/home/sturm/Downloads/adata_neutrophil_clusters.h5ad"
)

# %%
adata = sc.read_h5ad(
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/full_atlas_neutrophil_clusters.h5ad"
    # "/home/sturm/Downloads//full_atlas_neutrophil_clusters.h5ad"
)

# %%
adata_cpdb = sc.read_h5ad(
    "../../data/30_downstream_analyses/cell2cell_neutro/cpdb_h5ad/artifacts/adata_cpdb.h5ad"
)

# %%
patient_stratification = pd.read_csv(
    "../../data/30_downstream_analyses/stratify_patients/artifacts/patient_stratification.csv"
)

# %%
neutro_genesets = pd.read_excel(
    "../../tables/gene_annotations/neutro_phenotype_genesets.xlsx"
)
mask_valid_genes = neutro_genesets["gene_symbol"].isin(adata_n.var_names)
print(
    f"filtered out {np.sum(~mask_valid_genes)} genes because they are not in adata.var_names"
)
neutro_genesets = neutro_genesets.loc[mask_valid_genes].drop_duplicates()
neutro_genesets = {
    name: neutro_genesets.loc[lambda x: x["type"] == name, "gene_symbol"].tolist()
    for name in neutro_genesets["type"].unique()
}

# %%
neutro_recruitment_genes = pd.read_excel(
    "../../tables/gene_annotations/neutro_recruitment_chemokines.xlsx"
)

# %%
# statsmodels interprets "NAN" as np.nan, therefore let's keep "NAN-" and relabel later...
adata_n.obs["cell_type_tan_nan"] = [x[:4] for x in adata_n.obs["cell_type"]]

# %%
adata_n.obs["cell_type_tan_nan_label"] = [x[:3] for x in adata_n.obs["cell_type"]]

# %% [markdown]
# # UMAPs by covariate

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(
        adata_n,
        color="cell_type_tan_nan_label",
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
# sh.colors.set_scale_anndata(adata_n, "condition")
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(
        adata_n,
        color="condition",
        legend_fontoutline=2,
        frameon=False,
        size=20,
        add_outline=False,
        outline_width=[0.05, 0],
    )

# %%
with plt.rc_context({"figure.dpi": 150}):
    sc.pl.umap(
        adata_n[adata_n.obs["origin"] == "tumor_primary", :],
        color="condition",
        legend_fontoutline=2,
        frameon=False,
        size=20,
        title="condition: only cells from primary tumor samples",
    )

# %%
# sh.colors.set_scale_anndata(adata_n, "origin")
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
# # Quantify covariates by cluster

# %%
patients_with_neutros = (
    adata_n.obs.groupby("patient").size().loc[lambda x: x > 10].index
)

# %%
frac_by_patient = (
    adata_n.obs.groupby(["patient", "condition", "origin_biopsy"], observed=True)
    .apply(lambda x: x["cell_type_tan_nan_label"].value_counts(normalize=True))
    .unstack()
    .melt(ignore_index=False, var_name="cell_type", value_name="fraction")
    .reset_index()
    .loc[lambda x: x["patient"].isin(patients_with_neutros)]
)

# %%
frac_by_patient

# %%
alt.Chart(
    frac_by_patient.loc[lambda x: ~x["origin_biopsy"].isin(["nan", "tumor_metastasis"])]
).mark_boxplot().encode(x="cell_type", y="fraction", color="cell_type").facet(
    column="origin_biopsy"
)

# %%
alt.Chart(
    frac_by_patient.loc[lambda x: ~x["origin_biopsy"].isin(["nan", "tumor_metastasis"])]
).mark_boxplot().encode(x="cell_type", y="fraction", color="cell_type").facet(
    column="condition"
)

# %% [markdown]
# # Find marker genes for Neutrophil clusters

# %%
pb_tan_nan = sh.pseudobulk.pseudobulk(
    adata_n,
    groupby=[
        "cell_type_tan_nan_label",
        "patient",
        "condition",
        "dataset",
        "study",
        "tumor_stage",
    ],
)

# %%
sc.pp.normalize_total(pb_tan_nan, target_sum=1e6)
sc.pp.log1p(pb_tan_nan, base=2)

# %%
for ct in tqdm(pb_tan_nan.obs["cell_type_tan_nan_label"].unique()):
    sh.signatures.fold_change(
        pb_tan_nan,
        obs_col="cell_type_tan_nan_label",
        positive_class=ct,
        key_added=f"{ct}_fc",
    )
    sh.signatures.specific_fold_change(
        pb_tan_nan,
        obs_col="cell_type_tan_nan_label",
        positive_class=ct,
        key_added=f"{ct}_sfc",
    )
    sh.signatures.roc_auc(
        pb_tan_nan,
        obs_col="cell_type_tan_nan_label",
        positive_class=ct,
        key_added=f"{ct}_auroc",
    )

# %%
pb_n = sh.pseudobulk.pseudobulk(
    adata_n, groupby=["cell_type", "patient", "condition", "dataset", "tumor_stage"]
)

# %%
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n, base=2)

# %%
for ct in tqdm(pb_n.obs["cell_type"].unique()):
    sh.signatures.fold_change(
        pb_n, obs_col="cell_type", positive_class=ct, key_added=f"{ct}_fc"
    )
    sh.signatures.specific_fold_change(
        pb_n, obs_col="cell_type", positive_class=ct, key_added=f"{ct}_sfc"
    )
    sh.signatures.roc_auc(
        pb_n, obs_col="cell_type", positive_class=ct, key_added=f"{ct}_auroc"
    )


# %%
def metric_strip(
    pb, markers, *, metric_key="auroc", metric_title="AUROC", domain=[0.5, 1], top=10
):
    metric = []
    genes = []
    cts = []
    for ct, tmp_markers in markers.items():
        for gene in tmp_markers[:top]:
            metric.append(pb.var.loc[gene, f"{ct}_{metric_key}"])
            genes.append(gene)
            cts.append(ct)

    tmp_df = pd.DataFrame().assign(marker=genes, cell_type=cts, metric=metric)

    return (
        alt.Chart(tmp_df)
        .mark_rect()
        .encode(
            x=alt.X("marker", sort=None),
            color=alt.Color(
                "metric",
                scale=alt.Scale(scheme="viridis", domain=domain),
                title=metric_title,
            ),
        )
    )


def plot_markers(pb, groupby, markers, *, top=10):
    sc.pl.matrixplot(
        pb,
        groupby=groupby,
        var_names={k: v[:top] for k, v in markers.items()},
        cmap="bwr",
    )


# %%
markers_tan_nan = {
    ct: pb_tan_nan.var.loc[lambda x: (x[f"{ct}_auroc"] >= 0.75) & (x[f"{ct}_fc"] > 1.5)]
    .sort_values(f"{ct}_auroc", ascending=False)
    .index.tolist()
    for ct in sorted(pb_tan_nan.obs["cell_type_tan_nan_label"].unique())
}

# %%
markers = {
    ct: pb_n.var.loc[
        lambda x: (x[f"{ct}_auroc"] >= 0.75)
        & (x[f"{ct}_fc"] > 1.5)
        & (x[f"{ct}_sfc"] > 1)
    ]
    .sort_values(f"{ct}_auroc", ascending=False)
    .index.tolist()
    for ct in sorted(pb_n.obs["cell_type"].unique())
}

# %%
plot_markers(pb_tan_nan, "cell_type_tan_nan_label", markers_tan_nan, top=10)

# %%
plot_markers(pb_n, "cell_type", markers_tan_nan, top=10)

# %%
metric_strip(pb_tan_nan, markers_tan_nan)

# %%
plot_markers(pb_n, "cell_type", markers, top=5)

# %%
metric_strip(pb_n, markers, top=5)

# %% [markdown]
# ## Selected Top markers

# %%
selected_top_markers = {
    "NAN": ["CST7", "TMX4"],
    "TAN": ["CCL4", "CCL3"],
    "NAN-1": ["PADI4", "MMP9"],
    "NAN-2": ["CXCR2", "TNFSF13B"],
    "TAN-1": ["IL1A", "RHOH"],
    "TAN-2": ["CD74", "IFIT1"],
    "TAN-3": ["PLIN2", "PLPP3"],
    "TAN-4": ["RPL27", "RPS12"],
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

# %%
marker_res_tan_nan = (
    pd.concat(
        [
            pb_tan_nan.var.loc[lambda x: x["TAN_auroc"] > 0.5]
            .assign(group="TAN")
            .rename(columns={"TAN_fc": "fold-change", "TAN_auroc": "AUROC"}),
            pb_tan_nan.var.loc[lambda x: x["NAN_auroc"] > 0.5]
            .assign(group="NAN")
            .rename(columns={"NAN_fc": "fold-change", "NAN_auroc": "AUROC"}),
        ]
    )
    .loc[:, ["group", "fold-change", "AUROC"]]
    .sort_values("AUROC", ascending=False)
)
marker_res_tan_nan.to_csv(f"{artifact_dir}/tan_vs_nan_markers.csv")

# %%
pd.concat(
    [
        pb_n.var.loc[lambda x: x[f"{ct}_auroc"] > 0.5]
        .assign(group=ct)
        .rename(
            columns={
                f"{ct}_fc": "fold-change",
                f"{ct}_sfc": "specific-fold-change",
                f"{ct}_auroc": "AUROC",
            }
        )
        for ct in pb_n.obs["cell_type"].unique()
    ]
).loc[:, ["group", "fold-change", "specific-fold-change", "AUROC"]].sort_values(
    "AUROC", ascending=False
).to_csv(
    f"{artifact_dir}/neutro_cluster_markers.csv"
)

# %% [markdown]
# # Pairplot of candidate genes

# %%
pb_tan_nan.obs["cell_type_tan_nan_label"] = pd.Categorical(pb_tan_nan.obs["cell_type_tan_nan_label"] , categories=["NAN", "TAN"])

# %%
PROPS = {
    "boxprops": {"facecolor": "none", "edgecolor": "black"},
    "medianprops": {"color": "black"},
    "whiskerprops": {"color": "black"},
    "capprops": {"color": "black"},
}

genes_of_interest = [
        "CXCR4",
        "OLR1",
        "CD83",
        "VEGFA",
        "CXCR1",
        "CXCR2",
        "PTGS2",
        "SELL",
        "CSF3R",
    ]
fig = sh.pairwise.plot_paired(
    pb_tan_nan,
    "cell_type_tan_nan_label",
    paired_by="patient",
    var_names=genes_of_interest,
    hue="study",
    ylabel="log2(CPM+1)",
    size=6,
    n_cols=9,
    panel_size=(1.5,4),
    show=False,
    return_fig=True,
    pvalue_template=lambda x: "FDR<0.01" if x < 0.01 else f"FDR={x:.2f}",
    adjust_fdr=True,
    boxplot_properties=PROPS,
)
for i, ax in enumerate(fig.axes):
    ax.set_ylim(-0.5, np.max(pb_tan_nan[:, genes_of_interest].X) + 0.5)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if i > 0: 
        ax.yaxis.set_ticklabels([])
        ax.set_ylabel(None)

# %% [markdown] tags=[]
# # UMAP and Heatmaps of selected gene sets

# %% tags=[]
for gene_set, genes in neutro_genesets.items():
    display_html(f"<h2>Gene set of interest: {gene_set}</h2>", raw=True)
    sc.settings.set_figure_params(figsize=(2, 2))
    sc.pl.umap(
        adata_n,
        color=genes,
        cmap="inferno",
        size=10,
        ncols=10,
        frameon=False,
    )
    sc.pl.matrixplot(
        pb_n, var_names=genes, groupby=["cell_type"], cmap="bwr", title=gene_set
    )
    sc.pl.dotplot(adata_n, var_names=genes, groupby=["cell_type"], title=gene_set)
sc.settings.set_figure_params(figsize=(5, 5))

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
ct_fractions = ct_fractions.merge(
    patient_stratification, on=["patient", "study", "dataset"], how="left"
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

# %% tags=[]
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
ax.legend(bbox_to_anchor=(1.9, 1.05))
ax.set_title(
    "Neutrophil fraction in LUSC vs LUAD\np={:.4f}, linear model".format(
        res.pvalues["C(condition)[T.LUSC]"]
    )
)
plt.show()

# %% [markdown]
# ### Any differences between patient strata? 

# %%
neutro_subset2 = neutro_fractions.loc[
    lambda x: (
        x["dataset"].isin(datasets_with_neutros) & ~x["immune_infiltration"].isnull()
    ),
    :,
].sort_values("fraction", ascending=False)

# %%
neutro_subset2

# %%
mod = smf.ols(
    "fraction ~ C(immune_infiltration, Treatment('desert')) + dataset + tumor_type_annotated",
    data=neutro_subset2,
)
res = mod.fit()

# %%
res.summary()

# %%
fig, ax = plt.subplots()

neutro_subset2["immune_infiltration"] = neutro_subset2["immune_infiltration"].astype(
    str
)
neutro_subset2["dataset"] = neutro_subset2["dataset"].astype(str)

sns.stripplot(
    data=neutro_subset2,
    x="immune_infiltration",
    y="fraction",
    hue="study",
    palette=sh.colors.COLORS.study,
    ax=ax,
    size=7,
    linewidth=2,
)
ax.legend(bbox_to_anchor=(1.9, 1.05))
sns.boxplot(
    data=neutro_subset2,
    x="immune_infiltration",
    y="fraction",
    ax=ax,
    width=0.5,
    fliersize=0,
)
plt.show()

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
# # Neutrophil signatures
#
# Find marker genes that are specific for Neutrophils and specific for Neutrophil subclusters. 
# The purpose of the signatures is to use them for scoring bulk RNA-seq samples. 

# %%
neutro_sigs = {}

# %% [markdown]
# ### Compute metrics
#
# fold change, specific fold change and AUROC for each gene. This approach is defined in the MCP counter paper. 

# %%
sh.signatures.fold_change(
    pb_primary, obs_col="cell_type_major", positive_class="Neutrophils"
)

# %%
sh.signatures.specific_fold_change(
    pb_primary, obs_col="cell_type_major", positive_class="Neutrophils"
)

# %%
sh.signatures.roc_auc(
    pb_primary, obs_col="cell_type_major", positive_class="Neutrophils"
)

# %% [markdown]
# ## Neutrophil signatures
# Genes that are specific for the Neutrophil cluster as a whole

# %%
# Same principle as MCP counter, but even more specific!
# This is a subset of neutro_sig
tmp_df = (
    pb_primary.var.query("roc_auc >=0.97")
    .query("fold_change > 3")
    .query("specific_fold_change > 2")
)
neutro_sigs["neutro_sig"] = tmp_df.index.tolist()
tmp_df

# %%
# These are the original filtering criteria used by MCP counter
tmp_df = (
    pb_primary.var.query("roc_auc >=0.97")
    .query("fold_change >=2")
    .query("specific_fold_change >= 1.5")
)
neutro_sigs["neutro_sig2"] = tmp_df.index.tolist()
tmp_df

# %%
# Relaxed criteria to increase the number of genes, such that we can build TAN/NAN signatures.
tmp_df = (
    pb_primary.var.query("roc_auc >=0.90")
    .query("fold_change >=2")
    .query("specific_fold_change >= 1.5")
)
neutro_sigs["neutro_sig3"] = tmp_df.index.tolist()
tmp_df

# %% [markdown]
# ## TAN/NAN signatures
# Signature genes that are specific to Neutrophils *and* TANs or NANs. 
# We don't find genes that are specific enough to detect the TAN/NAN subclusters from bulk RNA-seq data. 

# %%
signature_genes_tan_nan = marker_res_tan_nan.loc[
    itertools.chain.from_iterable(markers_tan_nan.values())
]

# %% [markdown]
# ### Only nan genes in specificity-level 1 and 2

# %%
tmp_df = signature_genes_tan_nan.loc[lambda x: x.index.isin(neutro_sigs["neutro_sig"])]
tmp_df

# %%
tmp_df = signature_genes_tan_nan.loc[lambda x: x.index.isin(neutro_sigs["neutro_sig2"])]
neutro_sigs["nan_sig"] = tmp_df.loc[lambda x: x["group"] == "NAN"].index.tolist()
# neutro_sigs["tan_sig"] = tmp_df.loc[lambda x: x["group"] == "TAN"].index.tolist()
tmp_df

# %%
tmp_df = signature_genes_tan_nan.loc[
    lambda x: x.index.isin(neutro_sigs["neutro_sig3"]) & (x["group"] == "TAN")
]
neutro_sigs["tan_sig"] = tmp_df.index.tolist()
tmp_df

# %% [markdown]
# ## Evaluate signatures

# %%
{sig: len(genes) for sig, genes in neutro_sigs.items()}

# %% [markdown]
# ### Within Neutrophil cluster

# %%
for sig, genes in neutro_sigs.items():
    sc.tl.score_genes(pb_n, genes, score_name=sig)
    sc.tl.score_genes(adata_n, genes, score_name=sig)

# %%
sig_keys = list(neutro_sigs.keys()) + ["tan_nan_sig"]

# %%
adata_n.obs["tan_nan_sig"] = (adata_n.obs["tan_sig"] + adata_n.obs["nan_sig"]) / 2
pb_n.obs["tan_nan_sig"] = (pb_n.obs["tan_sig"] + pb_n.obs["nan_sig"]) / 2

# %%
sc.pl.matrixplot(pb_n, var_names=sig_keys, cmap="bwr", groupby="cell_type")

# %%
sc.pl.umap(
    adata_n,
    color=sig_keys,
    cmap="inferno",
    size=20,
    ncols=3,
    frameon=False,
)

# %% [markdown]
# ## Within all cell-types

# %%
for sig, genes in tqdm(neutro_sigs.items()):
    sc.tl.score_genes(pb_primary, genes, score_name=sig)
    sc.tl.score_genes(adata, genes, score_name=sig)

# %%
pb_primary.obs["tan_nan_sig"] = (
    pb_primary.obs["tan_sig"] + pb_primary.obs["nan_sig"]
) / 2

# %%
adata.obs["tan_nan_sig"] = (adata.obs["tan_sig"] + adata.obs["nan_sig"]) / 2

# %%
sc.pl.matrixplot(
    pb_primary,
    var_names=list(neutro_sigs.keys()) + ["tan_nan_sig"],
    cmap="bwr",
    groupby="cell_type_major",
)

# %%
sc.pl.umap(
    adata,
    color=list(neutro_sigs.keys()) + ["tan_nan_sig"],
    cmap="inferno",
    size=1,
    ncols=3,
    frameon=False,
)

# %% [markdown]
# ## Export signatures

# %%
with open(f"{artifact_dir}/neutro_sigs.csv", "w") as f:
    f.write("signature,gene_symbol\n")
    for sig, genes in neutro_sigs.items():
        for gene in genes:
            f.write(sig + "," + gene + "\n")

# %%

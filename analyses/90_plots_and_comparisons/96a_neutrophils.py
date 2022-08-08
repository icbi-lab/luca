# -*- coding: utf-8 -*-
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
from nxfvars import nxfvars
import sklearn.model_selection

# %%
alt.data_transformers.disable_max_rows()

# %%
sc.settings.set_figure_params(figsize=(5, 5))

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")
adata_n_path = nxfvars.get(
    "adata_n_path",
    # "../../data/30_downstream_analyses/neutrophils/subclustering/artifacts/adata_neutrophil_clusters.h5ad",
    "/home/sturm/Downloads/adata_neutrophil_clusters.h5ad",
)
adata_path = nxfvars.get(
    "adata_path",
    # "../../data/30_downstream_analyses/neutrophils/subclustering/artifacts/full_atlas_neutrophil_clusters.h5ad",
    "/home/sturm/Downloads/full_atlas_neutrophil_clusters.h5ad",
)
patient_stratification_path = nxfvars.get(
    "patient_stratification_path",
    "/home/sturm/Downloads/patient_stratification.csv"
    # "../../data/30_downstream_analyses/stratify_patients/stratification/artifacts/patient_stratification.csv",
)
neutro_geneset_path = nxfvars.get(
    "neutro_geneset_path",
    "../../tables/gene_annotations/neutro_phenotype_genesets.xlsx",
)
neutro_recruitment_geneset_path = nxfvars.get(
    "neutro_recruitment_geneset_path",
    "../../tables/gene_annotations/neutro_recruitment_chemokines.xlsx",
)
cpus = nxfvars.get("cpus", 8)

# %% [markdown]
# # Load data

# %%
adata_n = sc.read_h5ad(adata_n_path)

# %%
adata = sc.read_h5ad(adata_path)

# %%
patient_stratification = pd.read_csv(patient_stratification_path)

# %%
neutro_genesets = pd.read_excel(neutro_geneset_path)
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
neutro_recruitment_genes = pd.read_excel(neutro_recruitment_geneset_path)

# %%
ah = AnnotationHelper()

# %%
# statsmodels interprets "NAN" as np.nan, therefore let's keep "NAN-" and relabel later...
adata_n.obs["cell_type_tan_nan"] = [x[:4] for x in adata_n.obs["cell_type"]]

# %%
adata_n.obs["cell_type_tan_nan_label"] = [x[:3] for x in adata_n.obs["cell_type"]]

# %%
adata.obs["cell_type_major_tan_nan"] = adata.obs["cell_type_major"].astype(str)

# %%
adata.obs.loc[adata_n.obs_names, "cell_type_major_tan_nan"] = adata_n.obs[
    "cell_type_tan_nan"
]

# %% [markdown]
# ## Make Pseudobulk

# %%
pb_n = sh.pseudobulk.pseudobulk(
    adata_n, groupby=["cell_type", "patient", "condition", "dataset", "tumor_stage"]
)

# %%
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n, base=2)

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
pb_tan_nan.obs["cell_type_tan_nan_label"] = pd.Categorical(
    pb_tan_nan.obs["cell_type_tan_nan_label"], categories=["NAN", "TAN"]
)

# %%
sc.pp.normalize_total(pb_tan_nan, target_sum=1e6)
sc.pp.log1p(pb_tan_nan, base=2)

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true
# # UMAPs by covariate

# %%
with plt.rc_context({"figure.dpi": 150}):
    fig = sc.pl.umap(
        adata_n,
        color="cell_type_tan_nan_label",
        legend_loc="on data",
        legend_fontoutline=2,
        frameon=False,
        size=20,
        return_fig=True,
    )
    fig.savefig(
        f"{artifact_dir}/umap_neutro_tan_nan.pdf", dpi=1200, bbox_inches="tight"
    )

# %%
with plt.rc_context({"figure.dpi": 150}):
    fig = sc.pl.umap(
        adata_n,
        color="cell_type",
        legend_loc="on data",
        legend_fontoutline=2,
        frameon=False,
        size=20,
        return_fig=True,
    )
    fig.savefig(
        f"{artifact_dir}/umap_neutro_clusters.pdf", dpi=1200, bbox_inches="tight"
    )

# %%
# sh.colors.set_scale_anndata(adata_n, "condition")
with plt.rc_context({"figure.dpi": 150}):
    fig = sc.pl.umap(
        adata_n,
        color="condition",
        legend_fontoutline=2,
        frameon=False,
        size=20,
        add_outline=False,
        outline_width=[0.05, 0],
        return_fig=True,
    )
    fig.savefig(
        f"{artifact_dir}/umap_neutro_condition.pdf", dpi=1200, bbox_inches="tight"
    )

# %%
with plt.rc_context({"figure.dpi": 150}):
    fig = sc.pl.umap(
        adata_n[adata_n.obs["origin"] == "tumor_primary", :],
        color="condition",
        legend_fontoutline=2,
        frameon=False,
        size=20,
        title="condition: only cells from primary tumor samples",
        return_fig=True,
    )
    fig.savefig(
        f"{artifact_dir}/umap_neutro_condition_primary_tumor_only.pdf",
        dpi=1200,
        bbox_inches="tight",
    )

# %%
# sh.colors.set_scale_anndata(adata_n, "origin")
with plt.rc_context({"figure.dpi": 150}):
    fig = sc.pl.umap(
        adata_n[adata_n.obs["origin"] != "nan", :],
        color="origin",
        legend_fontoutline=2,
        frameon=False,
        size=20,
        return_fig=True,
    )
    fig.savefig(f"{artifact_dir}/umap_neutro_origin.pdf", dpi=1200, bbox_inches="tight")

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

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# # Clusters by patient and dataset

# %%
np.sum(adata_n.obs["dataset"].str.contains("UKIM-V"))

# %%
adata_n.shape[0]

# %%
np.sum(adata_n.obs["dataset"].str.contains("UKIM-V")) / adata_n.shape[0]

# %%
neutro_counts = (
    adata_n.obs.groupby(["cell_type", "dataset", "study", "patient"], observed=True)
    .size()
    .reset_index(name="n_cells")
    .groupby(["cell_type", "dataset"])
    .apply(lambda x: x.assign(n_cells_cell_type_dataset=lambda k: k["n_cells"].sum()))
    .groupby("patient")
    .apply(lambda x: x.assign(n_cells_patient=lambda k: k["n_cells"].sum()))
).query("n_cells_patient >= 10")

# %%
patient_cell_type_combs = (
    neutro_counts.loc[:, ["dataset", "study", "patient"]]
    .merge(neutro_counts.loc[:, ["cell_type"]], how="cross")
    .drop_duplicates()
)

# %%
tmp_df = neutro_counts.merge(
    patient_cell_type_combs,
    on=["dataset", "study", "patient", "cell_type"],
    how="outer",
)

# %%
tmp_df["n_cells"].fillna(0, inplace=True)

# %%
heatmp = (
    alt.Chart(tmp_df)
    .mark_rect()
    .encode(
        x="cell_type",
        y=alt.Y("patient", axis=None),
        color=alt.Color("n_cells", scale=alt.Scale(scheme="inferno", reverse=True)),
    )
)
txt = (
    alt.Chart(tmp_df)
    .mark_text()
    .encode(
        x="cell_type",
        y=alt.Y("patient", axis=None),
        text="n_cells",
        color=alt.condition(
            alt.datum.n_cells < 500, alt.value("black"), alt.value("white")
        ),
    )
)
studies = (
    alt.Chart(tmp_df)
    .mark_rect()
    .encode(
        y="patient",
        color=alt.Color("study", scale=sh.colors.altair_scale("study"), legend=None),
    )
)

# %%
(studies | (heatmp + txt).properties(width=300)).configure_concat(spacing=0)

# %% [markdown]
# ## more compact representation

# %%
tmp_df2 = (
    tmp_df.groupby(["study", "cell_type"], observed=True)
    .agg(n_cells=("n_cells", np.sum))
    .reset_index()
)

# %%
tmp_df2_study = (
    neutro_counts.groupby("study", observed=True)
    .agg(n_patients_ge_10_cells=("patient", lambda x: x.nunique()))
    .reset_index()
    .assign(x="#patients with ≥ 10 Neutrophils")
)

# %%
heatmp = (
    alt.Chart(tmp_df2)
    .mark_rect()
    .encode(
        x="cell_type",
        y=alt.Y("study", axis=None),
        color=alt.Color("n_cells", scale=alt.Scale(scheme="inferno", reverse=True)),
    )
)
txt = (
    alt.Chart(tmp_df2)
    .mark_text()
    .encode(
        x="cell_type",
        y=alt.Y("study", axis=None),
        text="n_cells",
        color=alt.condition(
            alt.datum.n_cells < 600, alt.value("black"), alt.value("white")
        ),
    )
)
studies = (
    alt.Chart(tmp_df2_study)
    .mark_rect()
    .encode(
        y=alt.Y("study"),
        x=alt.X("x", title=None),
        color=alt.Color("study", scale=sh.colors.altair_scale("study"), legend=None),
    )
)
studies_txt = (
    alt.Chart(tmp_df2_study)
    .mark_text()
    .encode(y="study", x=alt.X("x", title=None), text="n_patients_ge_10_cells")
)

# %%
studies + studies_txt

# %%
((studies + studies_txt) | (heatmp + txt).properties(width=300)).configure_concat(
    spacing=0
)

# %% [markdown]
# ## bar charts

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
patient_fracs.loc[lambda x: x["fraction"] > 0, "patient"].nunique()

# %%
ch = (
    alt.Chart(dataset_fracs)
    .mark_bar()
    .encode(
        x="fraction",
        y="cell_type",
        color=alt.Color("dataset", legend=alt.Legend(labelLimit=1e6)),
    )
)
ch.save(f"{artifact_dir}/neutro_clusters_by_dataset.svg")
ch.display()

# %%
ch = (
    alt.Chart(patient_fracs)
    .mark_bar()
    .encode(
        x=alt.X("fraction", scale=alt.Scale(domain=[0, 1])),
        y="cell_type",
        color=alt.Color("patient", legend=alt.Legend(labelLimit=1e6, symbolLimit=1000)),
    )
)
ch.save(f"{artifact_dir}/neutro_clusters_by_patient.svg")
ch.display()

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# # Quantify covariates by cluster

# %%
patients_with_neutros = (
    adata_n.obs.groupby("patient").size().loc[lambda x: x > 10].index
)

# %%
patients_with_neutros.nunique()

# %%
tan_nan_by_patient = (
    adata_n.obs.groupby(["patient", "condition", "origin_biopsy"], observed=True)
    .apply(lambda x: x["cell_type_tan_nan_label"].value_counts(normalize=True))
    .unstack()
    .melt(ignore_index=False, var_name="cell_type", value_name="fraction")
    .reset_index()
    .loc[lambda x: x["patient"].isin(patients_with_neutros)]
)

# %%
clusters_by_patient = (
    adata_n.obs.groupby(["patient", "condition", "origin_biopsy"], observed=True)
    .apply(lambda x: x["cell_type"].value_counts(normalize=True))
    .unstack()
    .melt(ignore_index=False, var_name="cell_type", value_name="fraction")
    .reset_index()
    .loc[lambda x: x["patient"].isin(patients_with_neutros)]
)

# %%
patient_meta = (
    adata.obs.loc[:, ["patient", "condition", "study"]]
    .drop_duplicates()
    .set_index("patient")
)

adata_tan_nan_by_patient = sc.AnnData(
    tan_nan_by_patient.pivot_table(
        values="fraction", index="patient", columns="cell_type"
    )
)
adata_tan_nan_by_patient.obs = adata_tan_nan_by_patient.obs.join(patient_meta)

adata_clusters_by_patient = sc.AnnData(
    clusters_by_patient.pivot_table(
        values="fraction", index="patient", columns="cell_type"
    )
)
adata_clusters_by_patient.obs = adata_clusters_by_patient.obs.join(patient_meta)

# %%
fig = sh.pairwise.plot_paired(
    adata_tan_nan_by_patient[
        adata_tan_nan_by_patient.obs["condition"].isin(["LUAD", "LUSC"]), :
    ],
    "condition",
    hue="study",
    n_cols=10,
    size=6,
    panel_size=(1.5, 4),
    pvalue_template=lambda x: f"p={x:.2f}",
    return_fig=True,
    show=False,
    ylabel="fraction",
    boxplot_properties={
        "boxprops": {"facecolor": "none", "edgecolor": "black"},
        "medianprops": {"color": "black"},
        "whiskerprops": {"color": "black"},
        "capprops": {"color": "black"},
    },
)
for i, ax in enumerate(fig.axes):
    ax.set_ylim(-0.05, 1.05)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if i > 0:
        ax.yaxis.set_ticklabels([])
        ax.set_ylabel(None)

# %%
fig = sh.pairwise.plot_paired(
    adata_clusters_by_patient[
        adata_clusters_by_patient.obs["condition"].isin(["LUAD", "LUSC"]), :
    ],
    "condition",
    hue="study",
    n_cols=10,
    size=6,
    panel_size=(1.5, 4),
    pvalue_template=lambda x: f"p={x:.2f}",
    return_fig=True,
    show=False,
    ylabel="fraction",
    boxplot_properties={
        "boxprops": {"facecolor": "none", "edgecolor": "black"},
        "medianprops": {"color": "black"},
        "whiskerprops": {"color": "black"},
        "capprops": {"color": "black"},
    },
)
for i, ax in enumerate(fig.axes):
    ax.set_ylim(-0.05, 1.05)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if i > 0:
        ax.yaxis.set_ticklabels([])
        ax.set_ylabel(None)

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[]
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

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true
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
fig.savefig(f"{artifact_dir}/neutro_fraction_luad_lusc.pdf", bbox_inches="tight")
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

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
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
fig.savefig(f"{artifact_dir}/vegfa_fractions.pdf", bbox_inches="tight")

# %% [markdown]
# # Pairplot of candidate genes

# %%
PROPS = {
    "boxprops": {"facecolor": "none", "edgecolor": "black"},
    "medianprops": {"color": "black"},
    "whiskerprops": {"color": "black"},
    "capprops": {"color": "black"},
}

genes_of_interest = {
    "TAN": ["OLR1", "VEGFA", "CD83", "ICAM1", "CXCR4"],
    "NAN": ["CXCR1", "CXCR2", "PTGS2", "SELL", "CSF3R", "FCGR3B"],
}
fig = sh.pairwise.plot_paired(
    pb_tan_nan,
    "cell_type_tan_nan_label",
    paired_by="patient",
    var_names=list(itertools.chain.from_iterable(genes_of_interest.values())),
    hue="study",
    ylabel="log2(CPM+1)",
    size=6,
    n_cols=12,
    panel_size=(1.5, 4),
    show=False,
    return_fig=True,
    pvalue_template=lambda x: "FDR<0.01" if x < 0.01 else f"FDR={x:.2f}",
    adjust_fdr=True,
    boxplot_properties=PROPS,
)
for i, ax in enumerate(fig.axes):
    ax.set_ylim(
        -0.5,
        np.max(
            pb_tan_nan[
                :, list(itertools.chain.from_iterable(genes_of_interest.values()))
            ].X
        )
        + 0.5,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if i > 0:
        ax.yaxis.set_ticklabels([])
        ax.set_ylabel(None)

fig.savefig(f"{artifact_dir}/pairplot_candidate_genes.pdf", bbox_inches="tight")

# %%
fig = plot_markers(pb_n, "cell_type", genes_of_interest, top=10, return_fig=True)
fig.savefig(f"{artifact_dir}/matrixplot_candidate_genes.pdf", bbox_inches="tight")

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# # Transcription factors

# %%
adata_n_dorothea = sh.compare_groups.compute_scores.run_dorothea(adata_n)

# %%
tf_res_tan_nan = sh.compare_groups.lm.lm_test_all(
    {"dorothea": adata_n_dorothea},
    groupby=["patient", "dataset"],
    column_to_test="cell_type_tan_nan",
    contrasts="Treatment('TAN-')",
    lm_covariate_str="+patient",
)

# %%
ch = sh.compare_groups.pl.plot_lm_result_altair(
    tf_res_tan_nan, cluster=True, p_cutoff=0.01
)
ch.save(f"{artifact_dir}/dorothea_tan_nan_heatmap_fdr_0.01.svg")
ch.display()

# %%
ch = sh.pairwise.plot_paired_fc(
    sh.pseudobulk.pseudobulk(
        adata_n_dorothea, groupby=["patient", "cell_type_tan_nan"], aggr_fun=np.mean
    ),
    var_names=tf_res_tan_nan.loc[
        lambda x: (x["fdr"] < 0.1) & (abs(x["coef"]) > 0.2), "variable"
    ].tolist(),
    metric="diff",
    groupby="cell_type_tan_nan",
    paired_by="patient",
)
ch.save(f"{artifact_dir}/dorothea_tan_nan_barchart.svg")
ch.display()

# %% [markdown]
# # Neutrophil subcluster signatures

# %%
adata_n_train, adata_n_test = sh.signatures.train_test_split(
    adata_n, replicate_col="patient"
)

# %% [markdown]
# ## TAN vs NAN

# %%
deseq_tan_nan = pd.read_csv(
    "../../data/30_downstream_analyses/de_analysis/tan_nan/de_deseq2/tan_nan_markers_adata_neutrophil_clusters_neutrophils_DESeq2_result.tsv",
    sep="\t",
)

# %%
deseq_tan_nan

# %%
markers_tan_nan = {
    "TAN": deseq_tan_nan.loc[
        lambda x: (x["log2FoldChange"] > 1) & (x["baseMean"] > 5), "gene_id"
    ].tolist(),
    "NAN": deseq_tan_nan.loc[
        lambda x: (x["log2FoldChange"] < -1) & (x["baseMean"] > 5), "gene_id"
    ].tolist(),
}

# %%
fig = sh.signatures.plot_markers(
    pb_tan_nan, "cell_type_tan_nan_label", markers_tan_nan, top=10, return_fig=True
)
fig.savefig(f"{artifact_dir}/matrixplot_tan_nan_top10.pdf", bbox_inches="tight")

# %%
fig = sh.signatures.plot_markers(
    pb_n, "cell_type", markers_tan_nan, top=10, return_fig=True, standard_scale="var"
)
fig.savefig(
    f"{artifact_dir}/matrixplot_tan_nan_top10_all_clusters.pdf", bbox_inches="tight"
)

# %%
pb_tan_nan.var["TAN_log_fdr"] = pb_tan_nan.var["NAN_log_fdr"] = -np.log10(
    deseq_tan_nan.set_index("gene_id")["padj"]
)

# %%
ch = sh.signatures.plot_metric_strip(
    pb_tan_nan,
    markers_tan_nan,
    metric_title="-log10(FDR)",
    metric_key="log_fdr",
    domain=[0, 30],
)
ch.save(f"{artifact_dir}/matrixplot_tan_nan_top10_log_fdr.svg")
ch.display()

# %%
fig = sh.signatures.plot_markers(
    pb_tan_nan, "cell_type_tan_nan_label", markers_tan_nan, top=30, return_fig=True
)
fig.savefig(f"{artifact_dir}/matrixplot_tan_nan_top30.pdf", bbox_inches="tight")

# %%
ch = sh.pairwise.plot_paired_fc(
    pb_tan_nan,
    "cell_type_tan_nan_label",
    var_names=markers_tan_nan["TAN"][:30] + markers_tan_nan["NAN"][:30],
    paired_by="patient",
    metric="diff",
    metric_name="log2 fold-change",
)
ch.save(f"{artifact_dir}/barchar_tan_nan_top30.svg")
ch.display()

# %% [markdown]
# ## Neutrophil clusters

# %%
deseq_neutro_clusters = pd.read_csv(
    "../../data/30_downstream_analyses/de_analysis/neutrophil_subclusters/de_deseq2/neutrophil_subclusters_adata_neutrophil_clusters_neutrophils_DESeq2_result.tsv",
    sep="\t",
)

# %%
deseq_neutro_clusters["log_fc_rank"] = deseq_neutro_clusters.groupby(["gene_id"])[
    "log2FoldChange"
].rank(ascending=False)

# %%
markers_neutro_clusters = {
    group: deseq_neutro_clusters.loc[
        lambda x: (x["log2FoldChange"] > 0)
        & (x["log_fc_rank"] == 1)
        & (x["baseMean"] > 5)
        & (x["comparison"] == group),
        "gene_id",
    ].tolist()
    for group in deseq_neutro_clusters["comparison"].unique()
}

# %%
fig = sh.signatures.plot_markers(
    pb_n,
    "cell_type",
    markers_neutro_clusters,
    top=5,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/matrixplot_neutro_clusters_top5.pdf", bbox_inches="tight")

# %% [markdown]
# ---

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
fig = plot_markers(
    pb_tan_nan, "cell_type_tan_nan_label", markers_tan_nan, top=10, return_fig=True
)
fig.savefig(f"{artifact_dir}/matrixplot_tan_nan_top10.pdf", bbox_inches="tight")

# %%
fig = plot_markers(
    pb_tan_nan, "cell_type_tan_nan_label", markers_tan_nan, top=30, return_fig=True
)
fig.savefig(f"{artifact_dir}/matrixplot_tan_nan_top30.pdf", bbox_inches="tight")

# %%
ch = sh.pairwise.plot_paired_fc(
    pb_tan_nan,
    "cell_type_tan_nan_label",
    var_names=markers_tan_nan["TAN"][:30] + markers_tan_nan["NAN"][:30],
    paired_by="patient",
    metric="diff",
    metric_name="log2 fold-change",
)
ch.save(f"{artifact_dir}/barchar_tan_nan_top30.svg")
ch.display()

# %%
fig = plot_markers(pb_n, "cell_type", markers_tan_nan, top=10, return_fig=True)
fig.savefig(
    f"{artifact_dir}/matrixplot_tan_nan_top10_all_clusters.pdf", bbox_inches="tight"
)

# %%
ch = metric_strip(pb_tan_nan, markers_tan_nan)
ch.save(f"{artifact_dir}/matrixplot_tan_nan_top10_auroc.svg")
ch.display()

# %%
np.max(pb_tan_nan.var.loc[:, pb_tan_nan.var.columns.str.contains("auroc")])

# %%
fig = plot_markers(pb_n, "cell_type", markers, top=5, return_fig=True)
fig.savefig(f"{artifact_dir}/matrixplot_neutro_clusters_top5.pdf", bbox_inches="tight")

# %%
ch = metric_strip(pb_n, markers, top=5)
ch.save(f"{artifact_dir}/matrixplot_neutro_clusters_top5_auroc.svg")
ch.display()

# %%
fig = plot_markers(pb_n, "cell_type", markers, top=10, return_fig=True)
fig.savefig(f"{artifact_dir}/matrixplot_neutro_clusters_top10.pdf", bbox_inches="tight")

# %%
ch = metric_strip(pb_n, markers, top=10)
ch.save(f"{artifact_dir}/matrixplot_neutro_clusters_top10_auroc.svg")
ch.display()

# %%
np.max(pb_n.var.loc[:, pb_n.var.columns.str.contains("auroc")])

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
selected_markers_of_interest = {
    "NAN-1": ["S100A12", "PADI4", "PROK2", "MMP9"],
    "TAN-1": ["CXCL8", "CXCL1", "CXCL2", "CD44"],
    "TAN-2": [
        "HLA-DRA",
        "CD74",
        "HLA-DPB1",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "XAF1",
        "ISG15",
    ],
    "TAN-3": ["PLIN2", "PLPP3", "MAP1LC3B", "PLAU"],
    "TAN-4": ["RPL10", "RPS2", "RPS18"],
}

# %%
fig = sc.pl.matrixplot(
    pb_n,
    var_names=selected_top_markers,
    groupby="cell_type",
    cmap="bwr",
    return_fig=True,
)
fig.savefig(
    f"{artifact_dir}/matrixplot_neutro_clusters_selected_markers.pdf",
    bbox_inches="tight",
)

# %%
fig = sc.pl.umap(
    adata_n,
    color=itertools.chain.from_iterable(selected_top_markers.values()),
    cmap="inferno",
    size=40,
    ncols=5,
    frameon=False,
    return_fig=True,
)
fig.savefig(
    f"{artifact_dir}/umap_neutro_clusters_selected_markers.pdf",
    bbox_inches="tight",
    dpi=1200,
)

# %%
fig = sc.pl.matrixplot(
    pb_n,
    var_names=selected_markers_of_interest,
    groupby="cell_type",
    cmap="bwr",
    return_fig=True,
)
fig.savefig(
    f"{artifact_dir}/matrixplot_neutro_clusters_selected_markers_of_interest.pdf",
    bbox_inches="tight",
)

# %%
fig = sc.pl.umap(
    adata_n,
    color=itertools.chain.from_iterable(selected_markers_of_interest.values()),
    cmap="inferno",
    size=40,
    ncols=5,
    frameon=False,
    return_fig=True,
)
fig.savefig(
    f"{artifact_dir}/umap_neutro_clusters_selected_markers_of_interest.pdf",
    bbox_inches="tight",
    dpi=300,
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
# # Neutrophil signatures
#
# Find marker genes that are specific for Neutrophils and specific for Neutrophil subclusters. 
# The purpose of the signatures is to use them for scoring bulk RNA-seq samples. 

# %%
neutro_sigs = {}

# %%
adata_primary = adata[adata.obs["origin"] == "tumor_primary", :]

# %%
adata_primary.var["marker"] = [
    "TAN"
    if g in markers_tan_nan["TAN"]
    else ("NAN" if g in markers_tan_nan["NAN"] else np.nan)
    for g in adata_primary.var.index
]

# %%
adata_primary_train, adata_primary_test = sh.signatures.train_test_split(
    adata_primary, replicate_col="patient"
)

# %%
results_neutro = sh.signatures.grid_search_cv(
    adata_primary_train,
    replicate_col="patient",
    label_col="cell_type_major",
    positive_class="Neutrophils",
    param_grid={
        "min_fc": list(np.arange(0.5, 3, 0.1)),
        "min_sfc": list(np.arange(0.5, 3, 0.1)),
        "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
    },
)

# %%
results_neutro.iloc[
    [
        0,
        1,
        2,
        3,
        4,
        50,
        100,
        200,
        300,
        400,
        500,
        600,
        700,
        800,
        900,
        1000,
        2000,
        2500,
        3000,
    ],
    :,
]

# %%
for i in [0, 1000, 2000, 3000]:
    tmp_mcpr = sh.signatures.refit_evaluate_plot(
        results_neutro,
        adata_train=adata_primary_train,
        adata_test=adata_primary_test,
        i=i,
    )
    neutro_sigs[f"sig_neutro_top_{i}"] = tmp_mcpr.signature_genes

# %%
for k, genes in list(neutro_sigs.items()):
    for ct in ["TAN", "NAN"]:
        if "sig_neutro" in k:
            neutro_sigs[k.replace("sig_neutro", f"sig_neutro_{ct.lower()}")] = [
                g for g in genes if g in markers_tan_nan[ct]
            ]


# %% [markdown]
# ## constrained to balanced

# %%
def balanced_constraint(var_df):
    var_df.sort_values("roc_auc", ascending=False, inplace=True)
    tan_genes = var_df.loc[lambda x: x["marker"] == "TAN"]
    nan_genes = var_df.loc[lambda x: x["marker"] == "NAN"]
    n = min(tan_genes.shape[0], nan_genes.shape[0])
    return pd.concat([tan_genes.iloc[:n, :], nan_genes.iloc[:n, :]])


# %%
results_neutro_balanced = sh.signatures.grid_search_cv(
    adata_primary_train,
    replicate_col="patient",
    label_col="cell_type_major",
    positive_class="Neutrophils",
    param_grid={
        "min_fc": list(np.arange(0.5, 3, 0.1)),
        "min_sfc": list(np.arange(0.5, 3, 0.1)),
        "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
        "constraint": [balanced_constraint]
    },
)

# %%
results_neutro_balanced.drop_duplicates(subset=["score_pearson", "n_genes"]).head(20)

 # %%
 tmp_mcpr = sh.signatures.refit_evaluate_plot(
    results_neutro_balanced,
    adata_train=adata_primary_train,
    adata_test=adata_primary_test,
    i=0,
)
neutro_sigs[f"sig_neutro_balanced"] = tmp_mcpr.signature_genes

# %% [markdown]
# ## constrained to TAN

# %%
results_neutro_constrained_tan = sh.signatures.grid_search_cv(
    adata_primary_train,
    replicate_col="patient",
    label_col="cell_type_major",
    positive_class="Neutrophils",
    param_grid={
        "min_fc": list(np.arange(0.5, 3, 0.1)),
        "min_sfc": list(np.arange(0.5, 3, 0.1)),
        "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
        "constraint": [lambda df: df.loc[lambda x: x["marker"] == "TAN"]]
    },
)

# %%
results_neutro_constrained_tan.iloc[[0, 50, 100, 200, 300, 500], :]

# %%
tmp_mcpr = sh.signatures.refit_evaluate_plot(
    results_neutro_constrained_tan,
    adata_train=adata_primary_train,
    adata_test=adata_primary_test,
    i=0,
)
neutro_sigs[f"sig_neutro_constrained_tan"] = tmp_mcpr.signature_genes


# %% [markdown]
# ## constrained to NAN

# %%
def constraint_nan(df):
    return df.loc[lambda x: x["marker"] == "NAN"]
    
results_neutro_constrained_nan = sh.signatures.grid_search_cv(
    adata_primary_train,
    replicate_col="patient",
    label_col="cell_type_major",
    positive_class="Neutrophils",
    param_grid={
        "min_fc": list(np.arange(0.5, 3, 0.1)),
        "min_sfc": list(np.arange(0.5, 3, 0.1)),
        "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
        "constraint": [constraint_nan]
    },
    n_jobs=cpus
)

# %%
results_neutro_constrained_nan.iloc[[0, 50, 100, 200, 300, 500], :]

# %%
tmp_mcpr = sh.signatures.refit_evaluate_plot(
    results_neutro_constrained_nan,
    adata_train=adata_primary_train,
    adata_test=adata_primary_test,
    i=0,
)
neutro_sigs[f"sig_neutro_constrained_nan"] = tmp_mcpr.signature_genes

# %% [markdown]
# # global TAN/NAN signatures

# %%
adata_primary_train.obs

# %%
results_tan = sh.signatures.grid_search_cv(
    adata_primary_train,
    replicate_col="patient",
    label_col="cell_type_major_tan_nan",
    positive_class="TAN-",
    param_grid={
        "min_fc": list(np.arange(0.5, 3, 0.1)),
        "min_sfc": list(np.arange(0.5, 3, 0.1)),
        "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
    },
)

# %%
results_tan

# %%
mcpr_tan = sh.signatures.refit_evaluate_plot(
    results_tan, adata_train=adata_primary_train, adata_test=adata_primary_test
)

# %%
neutro_sigs["sig_tan"] = mcpr_tan.signature_genes

# %% [markdown]
# ## NAN

# %%
results_nan = sh.signatures.grid_search_cv(
    adata_primary_train,
    replicate_col="patient",
    label_col="cell_type_major_tan_nan",
    positive_class="NAN-",
    param_grid={
        "min_fc": list(np.arange(0.5, 3, 0.1)),
        "min_sfc": list(np.arange(0.5, 3, 0.1)),
        "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
    },
)

# %%
results_nan

# %%
mcpr_nan = sh.signatures.refit_evaluate_plot(
    results_nan, adata_train=adata_primary_train, adata_test=adata_primary_test
)

# %%
neutro_sigs["sig_nan"] = mcpr_nan.signature_genes

# %% [markdown]
# ## Evaluate signatures

# %%
neutro_sigs = {k: v for k, v in neutro_sigs.items() if len(v)}

# %%
{sig: len(genes) for sig, genes in neutro_sigs.items()}

# %%
neutro_sigs["sig_trn"] = neutro_sigs["sig_tan"] + neutro_sigs["sig_nan"]

# %%
fig = sc.pl.matrixplot(
    pb_n,
    var_names=neutro_sigs["sig_neutro_top_0"],
    groupby="cell_type",
    cmap="bwr",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/neutro2_sig_genes_by_cluster.pdf", bbox_inches="tight")

# %% [markdown]
# ### Within Neutrophil cluster

# %%
for sig, genes in neutro_sigs.items():
    sc.tl.score_genes(pb_n, genes, score_name=sig)
    sc.tl.score_genes(adata_n, genes, score_name=sig)

# %%
# sig_keys = list(neutro_sigs.keys()) + ["tan_nan_sig"]

# %%
sig_keys = list(neutro_sigs.keys())

# %%
# adata_n.obs["tan_nan_sig"] = (adata_n.obs["tan_sig"] + adata_n.obs["nan_sig"]) / 2
# pb_n.obs["tan_nan_sig"] = (pb_n.obs["tan_sig"] + pb_n.obs["nan_sig"]) / 2

# %%
fig = sc.pl.matrixplot(
    pb_n, var_names=sig_keys, cmap="bwr", groupby="cell_type", return_fig=True
)
fig.savefig(
    f"{artifact_dir}/matrixplot_signature_scores_by_cluster.pdf", bbox_inches="tight"
)

# %%
fig = sc.pl.umap(
    adata_n,
    color=sig_keys,
    cmap="inferno",
    size=20,
    ncols=3,
    frameon=False,
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/umap_signature_scores.pdf", dpi=1200, bbox_inches="tight")

# %% [markdown]
# ## Within all cell-types

# %%
for sig, genes in tqdm(neutro_sigs.items()):
    sc.tl.score_genes(pb_primary, genes, score_name=sig)
    sc.tl.score_genes(adata, genes, score_name=sig)

# %%
sc.pl.matrixplot(
    pb_primary,
    var_names=list(neutro_sigs.keys()),  #  + ["tan_nan_sig"],
    cmap="bwr",
    groupby="cell_type_major",
)

# %%
fig = sc.pl.umap(
    adata,
    color=list(neutro_sigs.keys()),  #  + ["tan_nan_sig"],
    cmap="inferno",
    size=1,
    ncols=3,
    frameon=False,
    return_fig=True,
)
# fig.savefig(
#     f"{artifact_dir}/umap_atlas_signature_scores.pdf", dpi=1200, bbox_inches="tight"
# )

# %% [markdown]
# ## Export signatures

# %%
with open(f"{artifact_dir}/neutro_sigs.csv", "w") as f:
    f.write("signature,gene_symbol\n")
    for sig, genes in neutro_sigs.items():
        for gene in genes:
            f.write(sig + "," + gene + "\n")

# %%
assert False

# %% [markdown]
# ---

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
neutro_sigs = {}

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
adata_primary_train.obs["cell_type_major_tan_nan"].value_counts()

# %%
# results_tan_hierarch = sh.signatures.grid_search_cv(
#     adata_primary_train,
#     replicate_col="patient",
#     label_col="cell_type_major",
#     positive_class="Neutrophils",
#     label_col_fine="cell_type_major_tan_nan",
#     positive_class_fine="TAN-",
#     label_col_eval="cell_type_major_tan_nan",
#     positive_class_eval="TAN-",
#     model_class=sh.signatures.HierarchicalMCPSignatureRegressor,
#     param_grid={
#         "min_fc_coarse": list(np.arange(0.5, 3, 0.5)),
#         "min_sfc_coarse": list(np.arange(0.5, 3, 0.5)),
#         "min_auroc_coarse": [0.7, 0.8, 0.9, 0.95, 0.97],
#         "min_fc_fine": list(np.arange(0.5, 3, 0.5)),
#         "min_sfc_fine": [0],
#         "min_auroc_fine": [0.7, 0.8, 0.9, 0.95, 0.97],
#     },
# )

# %%
# pd.DataFrame(results_tan_hierarch).drop(columns="fold").groupby(
#     [
#         "min_fc_coarse",
#         "min_sfc_coarse",
#         "min_auroc_coarse",
#         "min_fc_fine",
#         "min_sfc_fine",
#         "min_auroc_fine",
#     ]
# ).agg(lambda x: x.mean(skipna=False)).sort_values(
#     "score_pearson", ascending=False
# ).head(
#     50
# )

# %%
# results_tan = sh.signatures.grid_search_cv(
#     adata_primary_train,
#     replicate_col="patient",
#     label_col="cell_type_major_tan_nan",
#     positive_class="TAN-",
#     param_grid={
#         "min_fc": list(np.arange(0.5, 3, 0.1)),
#         "min_sfc": list(np.arange(0.5, 3, 0.1)),
#         "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
#     },
# )

# %%
# pd.DataFrame(results_tan).drop(columns="fold").groupby(
#     ["min_fc", "min_sfc", "min_auroc"]
# ).agg(lambda x: x.mean(skipna=False)).sort_values(
#     "score_pearson", ascending=False
# ).head(
#     50
# )

# %%
# results_nan = sh.signatures.grid_search_cv(
#     adata_primary_train,
#     replicate_col="patient",
#     label_col="cell_type_major_tan_nan",
#     positive_class="NAN-",
#     param_grid={
#         "min_fc": list(np.arange(0.5, 3, 0.1)),
#         "min_sfc": list(np.arange(0.5, 3, 0.1)),
#         "min_auroc": [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.96, 0.97],
#     },
# )

# %%
# pd.DataFrame(results_nan).drop(columns="fold").groupby(
#     ["min_fc", "min_sfc", "min_auroc"]
# ).agg(lambda x: x.mean(skipna=False)).sort_values(
#     "score_pearson", ascending=False
# ).head(
#     50
# )

# %%
# mcpr_tan = sh.signatures.MCPSignatureRegressor(min_fc=2.8, min_sfc=0.5, min_auroc=0.97)
# sh.signatures.refit_and_evaluate(
#     mcpr_tan,
#     adata_primary_train,
#     adata_primary_test,
#     replicate_col="patient",
#     label_col="cell_type_major_tan_nan",
#     positive_class="TAN-",
# )

# %%
# mcpr_nan = sh.signatures.MCPSignatureRegressor(min_fc=2.7, min_sfc=1.7, min_auroc=0.85)
# sh.signatures.refit_and_evaluate(
#     mcpr_nan,
#     adata_primary_train,
#     adata_primary_test,
#     replicate_col="patient",
#     label_col="cell_type_major_tan_nan",
#     positive_class="NAN-",
# )

# %%
# neutro_sigs = {"NAN": mcpr_nan.signature_genes, "TAN": mcpr_tan.signature_genes}

# %%
# mcpr_nan.signature_genes

# %%
# mcpr_tan.signature_genes

# %%
# signature_genes_tan_nan = marker_res_tan_nan.loc[
#     itertools.chain.from_iterable(maarkers_tan_nan.values())
# ]

# %% [markdown]
# ## Only nan genes in specificity-level 1 and 2

# %%
# tmp_df = signature_genes_tan_nan.loc[lambda x: x.index.isin(neutro_sigs["neutro_sig"])]
# tmp_df

# %%
# tmp_df = signature_genes_tan_nan.loc[lambda x: x.index.isin(neutro_sigs["neutro_sig2"])]
# neutro_sigs["nan_sig"] = tmp_df.loc[lambda x: x["group"] == "NAN"].index.tolist()
# # neutro_sigs["tan_sig"] = tmp_df.loc[lambda x: x["group"] == "TAN"].index.tolist()
# tmp_df

# %%
# tmp_df = signature_genes_tan_nan.loc[
#     lambda x: x.index.isin(neutro_sigs["neutro_sig3"]) & (x["group"] == "TAN")
# ]
# neutro_sigs["tan_sig"] = tmp_df.index.tolist()
# tmp_df

# %%

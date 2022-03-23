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
)

# %%
adata = sc.read_h5ad(
    "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/full_atlas_neutrophil_clusters.h5ad"
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
print(f"filtered out {np.sum(~mask_valid_genes)} genes because they are not in adata.var_names")
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
pb_n = sh.pseudobulk.pseudobulk(
    adata_n, groupby=["cell_type", "patient", "condition", "dataset", "tumor_stage"]
)

# %%
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n, base=2)

# %%
pb_tan_nan = sh.pseudobulk.pseudobulk(
    adata_n,
    groupby=["cell_type_tan_nan", "patient", "condition", "dataset", "tumor_stage"],
)

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
        pb_tan_nan,
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
    .assign(group=lambda x: ["TAN" if _ > 0 else "NAN" for _ in x["log2_fc"]])
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
    var_names=top_markers_tan_nan["variable"],
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
sh.compare_groups.pl.plot_lm_result_altair(
    results["progeny"], p_cutoff=1, title="Progeny Neutrophils"
)

# %%
sh.compare_groups.pl.plot_lm_result_altair(
    results["dorothea"], title="Dorothea Neutrophils", p_cutoff=0.01
)

# %%
sh.compare_groups.pl.plot_lm_result_altair(
    results["cytosig"], title="Cytosig Neutrophils", p_cutoff=0.01
)

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
    sc.pl.matrixplot(pb_n, var_names=genes, groupby=["cell_type"], cmap="bwr", title=gene_set)
    sc.pl.dotplot(
        adata_n,
        var_names=genes,
        groupby=["cell_type"],
        title=gene_set
    )
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
ct_fractions = ct_fractions.merge(patient_stratification, on=["patient", "study", "dataset"], how="left")

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
mod = smf.ols("fraction ~ C(immune_infiltration, Treatment('desert')) + dataset + tumor_type_annotated", data=neutro_subset2)
res = mod.fit()

# %%
res.summary()

# %%
fig, ax = plt.subplots()

neutro_subset2["immune_infiltration"] = neutro_subset2["immune_infiltration"].astype(str)
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
    data=neutro_subset2, x="immune_infiltration", y="fraction", ax=ax, width=0.5, fliersize=0
)
plt.show()

# %% [markdown]
# # Cell2Cell communication

# %%
neutro_clusters = adata_n.obs["cell_type"].unique()

# %%
ad_cpdb_primary = adata_cpdb[
    (adata_cpdb.obs["origin"] == "tumor_primary"),
    (
        adata_cpdb.var["cluster_2"].isin(neutro_clusters)
        | adata_cpdb.var["cluster_1"].isin(neutro_clusters)
    ),
]
# keep only interactions that are >0 in at least 5 samples
ad_cpdb_primary = ad_cpdb_primary[:, np.sum(ad_cpdb_primary.X != 0, axis=0) > 5].copy()


# %%
def cpdb_for_cell_types(adata_cpdb, cell_types, direction="incoming"):
    """
    Regroup the CPDB anndata object to consider cell-types as samples (which allows
    comparisons between cell-types, rather than only patients.

    direction can be `incoming` for signalling towards the selected cell-types, or
    `outgoing` for signalling originating from the selected cell-types.
    """
    cpdb_molten = pd.DataFrame(
        adata_cpdb.X,
        index=adata_cpdb.obs_names,
        columns=adata_cpdb.var_names,
    ).melt(ignore_index=False)

    tmp_var = (
        adata_cpdb.var.loc[
            lambda x: x["cluster_1"].isin(cell_types)
            != x["cluster_2"].isin(cell_types)  # != is xor
        ]
        .assign(
            ligrec=lambda _: [
                ("incoming" if x["cluster_2"] in cell_types else "outgoing")
                for idx, x in _.iterrows()
            ]
        )
        .reset_index()
        .assign(
            idx_new=lambda x: [
                idx.replace("_" + (clus1 if ligrec == "outgoing" else clus2), "")
                for idx, ligrec, clus1, clus2 in zip(
                    x["idx"], x["ligrec"], x["cluster_1"], x["cluster_2"]
                )
            ]
        )
    )

    cpdb_molten2 = (
        cpdb_molten.reset_index()
        .set_index("idx")
        .rename(columns={"index": "patient"})
        .join(tmp_var.set_index("idx"), how="inner")
    )

    cpdb_molten_directed = cpdb_molten2.loc[lambda x: x["ligrec"] == direction]

    cluster_to_drop = "cluster_2" if direction == "incoming" else "cluster_1"
    cluster_to_keep = "cluster_1" if direction == "incoming" else "cluster_2"
    adata_directed = sc.AnnData(
        cpdb_molten_directed.reset_index()
        .assign(
            sample=lambda x: x["patient"].astype(str)
            + "_"
            + x[cluster_to_drop].astype(str),
        )
        .pivot(
            index=["patient", "sample", cluster_to_drop],
            columns="idx_new",
            values="value",
        )
        .fillna(0)
    )

    adata_directed.var = adata_directed.var.join(
        tmp_var.loc[
            lambda x: x["ligrec"] == direction,
            ["idx_new", "source", "target", cluster_to_keep],
        ]
        .set_index("idx_new")
        .drop_duplicates(),
        how="left",
    )
    adata_directed.obs = adata_directed.obs.reset_index().rename(
        columns={cluster_to_drop: "group"}
    )
    return adata_directed


# %%
adata_cpdb_neutro_incoming = cpdb_for_cell_types(
    ad_cpdb_primary, neutro_clusters, "incoming"
)

# %%
adata_cpdb_neutro_outgoing = cpdb_for_cell_types(
    ad_cpdb_primary, neutro_clusters, "outgoing"
)

# %%
res_incoming = sh.compare_groups.lm.test_lm(
    adata_cpdb_neutro_incoming, "~ C(group, Sum) + patient", "group"
)

# %%
res_outgoing = sh.compare_groups.lm.test_lm(
    adata_cpdb_neutro_outgoing, "~ C(group, Sum) + patient", "group"
)

# %%
res_incoming.loc[lambda x: x["variable"].str.contains("Tumor cells")].pipe(
    sh.util.fdr_correction
).pipe(sh.compare_groups.pl.plot_lm_result_altair, p_cutoff=0.01)

# %%
res_outgoing.loc[lambda x: x["variable"].str.contains("Tumor cells")].pipe(
    sh.util.fdr_correction
).pipe(sh.compare_groups.pl.plot_lm_result_altair, p_cutoff=0.01)

# %%
res_outgoing.loc[lambda x: x["variable"].str.contains("T cell CD8+")].pipe(
    sh.util.fdr_correction
).pipe(sh.compare_groups.pl.plot_lm_result_altair, p_cutoff=0.01)

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
    .query("specific_fold_change >= 1")
)
neutro_sigs["neutro_sig3"] = tmp_df.index.tolist()
tmp_df

# %% [markdown]
# ## TAN/NAN signatures
# Signature genes that are specific to Neutrophils *and* TANs or NANs. 
# We don't find genes that are specific enough to detect the TAN/NAN subclusters from bulk RNA-seq data. 

# %%
signature_genes_tan_nan = marker_res_tan_nan.loc[
    lambda x: (x["fdr"] < 0.01) & (abs(x["log2_fc"]) > 1)
]

# %% [markdown]
# ### Only nan genes in specificity-level 1 and 2

# %%
tmp_df = signature_genes_tan_nan.loc[
    lambda x: x["variable"].isin(neutro_sigs["neutro_sig"])
]
neutro_sigs["nan_sig"] = tmp_df["variable"].tolist()
tmp_df

# %%
tmp_df = signature_genes_tan_nan.loc[
    lambda x: x["variable"].isin(neutro_sigs["neutro_sig3"]) & (x["group"] == "TAN")
]
neutro_sigs["tan_sig"] = tmp_df["variable"].tolist()
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
sc.pl.matrixplot(
    pb_n, var_names=list(neutro_sigs.keys()), cmap="bwr", groupby="cell_type"
)

# %%
sc.pl.umap(
    adata_n,
    color=list(neutro_sigs.keys()),
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
sc.pl.matrixplot(
    pb_primary,
    var_names=list(neutro_sigs.keys()),
    cmap="bwr",
    groupby="cell_type_major",
)

# %%
sc.pl.umap(
    adata,
    color=list(neutro_sigs.keys()),
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

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
from hierarchical_bootstrapping.util import gini_index

# %%
import scanpy as sc
import pandas as pd
from nxfvars import nxfvars

sc.settings.set_figure_params(figsize=(5, 5))
from pathlib import Path
from scanpy_helpers.annotation import AnnotationHelper
import progeny
import dorothea
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
import hierarchical_bootstrapping as hb
import altair as alt

# %%
threadpool_limits(32)

# %%
ah = AnnotationHelper()

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
adata_epi = sc.read_h5ad(
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/33_cell_types_epi/artifacts/adata_epithelial.h5ad"
)

# %%
adata_tumor = sc.read_h5ad(
    "../../data/20_integrate_scrnaseq_data/annotate_datasets/33_cell_types_epi/artifacts/adata_tumor.h5ad"
)

# %%
sc.pl.umap(adata, color="cell_type")

# %%
scissor_res_files = Path(
    "../../data/20_integrate_scrnaseq_data/scissor/scissor_by_sample/artifacts/"
).glob("*_scissor.tsv")

# %%
scissor_ids = [pd.read_csv(x, sep="\t") for x in scissor_res_files]

# %%
scissor_obs = (
    pd.concat(scissor_ids)
    .set_index("cell_id")
    .rename(columns={"Scissor_select": "scissor"})
)

# %%
scissor_obs

# %%
adata.obs["scissor"] = scissor_obs["scissor"]

# %%
sc.settings.set_figure_params(figsize=(8, 8))

# %%
sc.pl.umap(adata, color=["scissor", "cell_type"], size=1)

# %%
sc.pl.umap(adata, color="VEGFA")

# %%
sc.pl.umap(adata, color=["scissor", "cell_type_coarse"], size=1)

# %%
sc.pl.umap(adata, color=["dataset", "condition", "origin"], size=1)

# %%
adata_epi.obs["scissor"] = adata.obs["scissor"]
adata_epi.obs["cell_type"] = adata.obs["cell_type"]
adata_tumor.obs["scissor"] = adata.obs["scissor"]

# %%
sc.pl.umap(adata_epi, color="cell_type")

# %%
sc.pl.umap(adata_epi, color=["dataset", "condition", "origin"], size=2, wspace=0.5)

# %%
sc.pl.umap(adata_epi, color=["scissor", "cell_type"], size=2, wspace=0.5)

# %%
sc.pl.umap(
    adata_epi[adata_epi.obs["origin"] == "tumor_primary", :], color="scissor", size=4
)

# %%
sc.pl.umap(adata_tumor, color=["scissor", "cell_type"], size=4, wspace=0.5)

# %%
sc.pl.umap(
    adata[adata.obs["origin"] == "tumor_primary", :],
    color=["scissor", "cell_type"],
    size=4,
)

# %%
adata_primary = adata[adata.obs["origin"] == "tumor_primary", :].copy()

# %%
adata_primary.obs["scissor"] = adata_primary.obs["scissor"].astype(str)

# %%
scissor_per_cell_type = (
    (
        adata_primary.obs.groupby(["cell_type"])["scissor"]
        .value_counts(normalize=True)
        .reset_index(name="frac")
    )
    .pivot_table(values="frac", columns="scissor", index="cell_type", fill_value=0)
    .reset_index()
)

# %%
scissor_per_cell_type

# %%
# TODO aggregate by patient

# %%
alt.Chart(scissor_per_cell_type).mark_bar().encode(
    x=alt.X("cell_type", sort="y"), y="better survival"
)

# %%
alt.Chart(scissor_per_cell_type).mark_bar(color="orange").encode(
    x=alt.X("cell_type", sort="-y"), y="worse survival"
)


# %%
def scissor_gini(tmp_adata):
    # tmp_adata = adata[adata.obs["cell_type"] == cell_type, :].copy()
    sc.pp.neighbors(tmp_adata, use_rep="X_scANVI")
    sc.tl.leiden(tmp_adata, resolution=0.5)
    fractions = (
        tmp_adata.obs.groupby("leiden")["scissor"]
        .value_counts(normalize=True)
        .reset_index(name="frac")
        .pivot_table(values="frac", columns="scissor", index="leiden", fill_value=0)
        .reset_index()
    )
    try:
        gini_better = gini_index(fractions["better survival"].values)
    except KeyError:
        gini_better = 0
    try:
        gini_worse = gini_index(fractions["worse survival"].values)
    except KeyError:
        gini_worse = 0
    return gini_better, gini_worse


# %%
from tqdm.contrib.concurrent import process_map
import itertools

# %%
adatas = []
cell_types = adata_primary.obs["cell_type"].unique()
for cell_type in cell_types:
    adatas.append(adata_primary[adata_primary.obs["cell_type"] == cell_type, :].copy())

# %%
gini_res = process_map(scissor_gini, adatas)

# %%
gini_better, gini_worse = zip(*gini_res)

# %%
gini_better = pd.Series(gini_better, index=cell_types)
gini_worse = pd.Series(gini_worse, index=cell_types)

# %%
scissor_per_cell_type.set_index("cell_type", inplace=True)

# %%
scissor_per_cell_type["gini_better"] = gini_better
scissor_per_cell_type["gini_worse"] = gini_worse
scissor_per_cell_type = scissor_per_cell_type.reset_index()

# %%
scissor_per_cell_type

# %%
order = scissor_per_cell_type.sort_values("better survival")[
    "cell_type"
].values.tolist()

# %%

# %%
alt.vconcat(
    alt.Chart(scissor_per_cell_type)
    .mark_bar()
    .encode(
        x=alt.X("cell_type", sort=order, axis=None),
        y="better survival",
        color=alt.Color("gini_better", scale=alt.Scale(scheme="magma", reverse=False)),
    )
    .properties(height=100),
    alt.Chart(
        scissor_per_cell_type.assign(
            **{"worse survival": lambda x: -x["worse survival"]}
        )
    )
    .mark_bar()
    .encode(
        x=alt.X("cell_type", sort=order),
        y="worse survival",
        color=alt.Color("gini_worse", scale=alt.Scale(scheme="magma", reverse=False)),
    )
    .properties(height=100),
    spacing=0,
)

# %% [markdown]
# ---

# %% [markdown]
# # Patient stratification

# %%
adata_primary = adata[
    (adata.obs["origin"] == "tumor_primary")
    # exclude datasets that only contain a single cell-type
    & ~adata.obs["dataset"].isin(["Guo_Zhang_2018_NSCLC", "Maier_Merad_2020_NSCLC"]),
    :,
]

# %%
cell_type_fractions = (
    adata_primary.obs.groupby(["patient", "dataset", "condition"], observed=True)[
        "cell_type"
    ]
    .value_counts(normalize=True)
    .reset_index(name="frac")
)

# %%
ad_ctf = sc.AnnData(
    X=cell_type_fractions.pivot(
        index=["patient", "dataset", "condition"], columns="level_3", values="frac"
    )
)

# %%
ad_ctf.obs = ad_ctf.obs.reset_index()

# %%
sc.pp.regress_out(ad_ctf, "dataset")

# %%
sc.tl.dendrogram(ad_ctf, groupby="patient", use_rep="X")

# %%
sc.pl.matrixplot(
    ad_ctf, var_names=ad_ctf.var_names, groupby="patient", dendrogram=True, swap_axes=True, cmap="bwr", vmin=-.5, vmax=.5
)

# %%
sc.pp.neighbors(ad_ctf, use_rep="X")

# %%
sc.tl.umap(ad_ctf)

# %%
sc.pl.umap(ad_ctf, color=ad_ctf.var_names,cmap="bwr", vmin=-.5, vmax=.5, add_outline=True)

# %%
sc.tl.leiden(ad_ctf)

# %%
sc.pl.umap(ad_ctf, color=["leiden", "dataset", "condition"])

# %%
sc.tl.rank_genes_groups(ad_ctf, groupby="leiden", method="wilcoxon")

# %%
sc.pl.rank_genes_groups_matrixplot(ad_ctf, values_to_plot="log10_pvals_adj", dendrogram=False, n_genes=4, vmax=5)

# %%
sc.pl.matrixplot(ad_ctf, var_names=ad_ctf.var_names, groupby="leiden", cmap="bwr", vmax=0.25, vmin=-0.25, swap_axes=True)

# %%
ah.annotate_cell_types(ad_ctf, {
    "T/B cell infiltration": [0, 5],
    "neutral": [1],
    "immune-excluded": [2,3],
    "myeloid infiltration": [4],
    "Club": [6]
}, key_added="patient_group")

# %%
sc.pl.matrixplot(ad_ctf, var_names=ad_ctf.var_names, groupby="patient_group", cmap="bwr", vmax=0.25, vmin=-0.25, swap_axes=False)

# %%
patient_map = {p: g for p, g in zip(ad_ctf.obs["patient"], ad_ctf.obs["patient_group"])}

# %%
adata.obs["patient_group"] = [patient_map.get(p, None) for p in adata.obs["patient"]]

# %%
sc.pl.umap(adata, color=["patient_group", "cell_type"])

# %% [markdown]
# ---

# %%
sc.settings.set_figure_params(figsize=(4, 4))

# %% tags=[]
regulons = dorothea.load_regulons(
    [
        "A",
        "B",
    ],  # Which levels of confidence to use (A most confident, E least confident)
    organism="Human",  # If working with mouse, set to Mouse
)

# %%
dorothea.run(
    adata,  # Data to use
    regulons,  # Dorothea network
    center=True,  # Center gene expression by mean per cell
    num_perm=0,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # TF with less than 5 targets will be ignored
)

# %%
adata_dorothea = dorothea.extract(adata)

# %%
model = progeny.load_model(
    organism="Human",  # If working with mouse, set to Mouse
    top=1000,  # For sc we recommend ~1k target genes since there are dropouts
)

# %%
progeny.run(
    adata,  # Data to use
    model,  # PROGENy network
    center=True,  # Center gene expression by mean per cell
    num_perm=0,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # Pathways with less than 5 targets will be ignored
)

# %%
adata_progeny = progeny.extract(adata)

# %%
cytosig_signature = pd.read_csv(
    "../../tables/cytosig_signature_matrix.tsv", sep="\t", index_col=0
)

# %%
progeny.run(
    adata,  # Data to use
    cytosig_signature,  # PROGENy network
    center=True,  # Center gene expression by mean per cell
    num_perm=0,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # Pathways with less than 5 targets will be ignored
    obsm_key="cytosig",
)

# %%
adata_cytosig = progeny.extract(adata, obsm_key="cytosig")

# %%
sc.pl.matrixplot(
    adata_cytosig,
    var_names=adata_cytosig.var_names,
    groupby="cell_type",
    cmap="bwr",
    vmin=-3,
    vmax=3,
    swap_axes=True,
)

# %%
sc.pl.matrixplot(
    adata_progeny,
    var_names=adata_progeny.var_names,
    groupby="cell_type",
    cmap="bwr",
    vmin=-3,
    vmax=3,
    swap_axes=True,
)

# %%
sc.pl.matrixplot(
    adata_dorothea,
    var_names=adata_dorothea.var_names,
    groupby="cell_type",
    cmap="bwr",
    vmin=-3,
    vmax=3,
    swap_axes=True,
)

# %%
import warnings
import scipy.stats


# %%
def paired_test(tmp_adata):
    var_names = tmp_adata.var_names
    groups = ["tumor", "normal"]
    tmp_adata.obs["status"] = [
        {"tumor_primary": "tumor", "normal": "normal", "normal_adjacent": "normal"}.get(
            origin, None
        )
        for origin in tmp_adata.obs["origin"]
    ]
    pseudobulk = hb.tl.pseudobulk(
        tmp_adata, groupby=["patient", "sample", "status"], aggr_fun=np.mean
    )
    paired_by = "patient"
    group_by = "status"
    df = pseudobulk.obs.loc[:, [paired_by, group_by]].join(
        pd.DataFrame(pseudobulk.X, index=pseudobulk.obs_names, columns=var_names)
    )
    # remove unpaired samples
    has_matching_samples = df.groupby(paired_by).apply(
        lambda x: sorted(x[group_by]) == sorted(groups)
    )
    has_matching_samples = has_matching_samples.index[has_matching_samples].values
    removed_samples = adata.obs[paired_by].nunique() - len(has_matching_samples)
    if removed_samples:
        warnings.warn(f"{removed_samples} unpaired samples removed")

    values1 = df.loc[
        df[paired_by].isin(has_matching_samples) & (df[group_by] == groups[0]),
        var_names,
    ]
    values2 = df.loc[
        df[paired_by].isin(has_matching_samples) & (df[group_by] == groups[1]),
        var_names,
    ]

    mean1 = np.mean(values1, axis=0)
    mean2 = np.mean(values2, axis=0)
    pvalues = []
    for col in var_names:
        _, p = scipy.stats.wilcoxon(values1[col].values, values2[col].values)
        pvalues.append(p)
    # _, pvalues = scipy.stats.ttest_rel(
    #     values1, values2
    # )

    return pd.DataFrame(
        index=var_names, data={"mean1": mean1, "mean2": mean2, "pvalue": pvalues}
    )


# %%
def run_progeny(adata):
    tmp_adata = adata.copy()
    progeny.run(
        tmp_adata,  # Data to use
        model,  # PROGENy network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # Pathways with less than 5 targets will be ignored
    )
    return progeny.extract(tmp_adata)


# %%
adatas_by_cell_type = []
cell_types = adata.obs["cell_type"].unique()
for cell_type in cell_types:
    adatas_by_cell_type.append(adata[adata.obs["cell_type"] == cell_type, :].copy())

# %% tags=[]
adatas_progeny = list(
    tqdm(map(run_progeny, adatas_by_cell_type), total=len(adatas_by_cell_type))
)

# %% tags=[]
cell_types

# %%
dfs = []
for ct, ad in zip(cell_types, adatas_progeny):
    try:
        dfs.append(paired_test(ad).assign(cell_type=ct))
    except ValueError:
        pass

# %%
pd.set_option("display.max_rows", None)

# %%
pd.concat(dfs).sort_values("pvalue")

# %%
pseudobulk_epi = 

# %%

# %%

# %%
from numba import njit

# %%
from hierarchical_bootstrapping.util import gini_index
import numpy as np
from tqdm import tqdm

# %%
progeny_epi = adata_cytosig[
    adata_cytosig.obs["cell_type_coarse"] == "Epithelial cell", :
].copy()

# %%
pb_tn = hb.tl.pseudobulk(
    progeny_epi, groupby=["dataset", "patient", "status"], aggr_fun=np.mean
)

# %%
hb.pl.paired_expression(
    pb_tn, group_by="status", paired_by="patient", var_names=None, show_legend=True
)

# %% [markdown]
# ---

# %%
# %%time
hb.tl.pseudobulk(progeny_epi, groupby=["leiden", "patient"], aggr_fun=np.mean)

# %%
progeny_epi_pb_leiden = hb.tl.pseudobulk(
    progeny_epi, groupby=["patient", "leiden"], aggr_fun=np.mean
)

# %%
progeny_epi_pb = hb.tl.pseudobulk(progeny_epi, groupby=["patient"], aggr_fun=np.mean)

# %%
progeny_epi_pb.X -= np.min(progeny_epi_pb.X)
progeny_epi_pb_leiden.X -= np.min(progeny_epi_pb_leiden.X)

# %%
res = []
patients = progeny_epi_pb_leiden.obs["patient"].unique()
for patient in tqdm(patients):
    tmp_x = progeny_epi_pb_leiden.X[progeny_epi_pb_leiden.obs["patient"] == patient, :]
    res.append(np.apply_along_axis(gini_index, 0, tmp_x))

# %%
gini_within = np.mean(np.vstack(res), axis=0)

# %%
gini_between = np.apply_along_axis(gini_index, 0, progeny_epi_pb.X)

# %%
df_to_plot = (
    pd.DataFrame(gini_within.T, index=progeny_epi.var_names, columns=["gini_within"])
    .join(
        pd.DataFrame(
            gini_between, index=progeny_epi.var_names, columns=["gini_between"]
        )
    )
    .reset_index()
)

# %%
df_to_plot

# %%
import altair as alt

# %%
alt.Chart(df_to_plot).mark_point().encode(
    y="gini_within", x="gini_between", color=df_to_plot.columns[0], tooltip="index"
)

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    ah.plot_umap(
        adata,
        filter_cell_type=["Goblet"],
        cmap="inferno",
    )

# %%
with plt.rc_context({"figure.figsize": (4, 4)}):
    ah.plot_umap(
        adata_epi,
        filter_cell_type=[
            "Ciliated",
            "Alevolar",
            "Basal",
            "Club",
            "Dividing",
            "Goblet",
            "Ionocyte",
            "Mesothelial",
            "Suprabasal",
        ],
        cmap="inferno",
    )

# %%
adata_goblet = adata_epi[adata_epi.obs["cell_type"] == "Ciliated", :]

# %%
adata_goblet

# %%
sc.pl.umap(adata_goblet, color=["dataset", "scissor"])

# %%
adata_m = adata[adata.obs["cell_type"].isin(["DC mature/cDC 1", "cDC2", "Monocyte"]), :]

# %%

# %%
sc.pl.umap(adata_m, color=["dataset", "cell_type", "scissor"])

# %%
ah.plot_umap(
    adata_m, filter_cell_type=["Macro", "Mono", "DC", "Div"], cmap="inferno", size=2
)

# %%

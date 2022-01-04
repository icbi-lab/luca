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
import re
import statsmodels.stats
import warnings
import scipy.stats
from tqdm.auto import tqdm
import numpy as np
import statsmodels.formula.api as smf
from tqdm.contrib.concurrent import process_map
import itertools

# %%
threadpool_limits(32)

# %%
ah = AnnotationHelper()

# %%
path_adata = nxfvars.get(
    "adata_in",
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)

# %%
adata = sc.read_h5ad(path_adata)

# %%
adata_epi = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/33_cell_types_epi/artifacts/adata_epithelial.h5ad"
)

# %%
adata_tumor = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/33_cell_types_epi/artifacts/adata_tumor.h5ad"
)

# %%
sc.pl.umap(adata, color="cell_type")

# %%
scissor_res_files = {
    id: Path(
        "../../data/30_downstream_analyses/scissor/scissor_by_sample/artifacts/"
    ).glob(f"*_scissor_{id}.tsv")
    for id in [
        # TODO something is wrong with e.g. lambrechts 6149v2_3: the cell ids get lost, but only for survival! 
        # "survival",
        "tumor_stage",
        "kras_mutation",
        "braf_mutation",
        "egfr_mutation",
    ]
}

# %%
scissor_ids = {
    id: [pd.read_csv(x, sep="\t") for x in tmp_files]
    for id, tmp_files in scissor_res_files.items()
}

# %%
scissor_obs = {
    f"scissor_{id}": (
        pd.concat(tmp_ids)
        .set_index("cell_id")
        .rename(columns={"Scissor_select": "scissor"})
    )
    for id, tmp_ids in scissor_ids.items()
}

# %%
for colname, series in scissor_obs.items():
    print(colname)
    adata.obs[colname] = series

# %%
sc.settings.set_figure_params(figsize=(8, 8))

# %%
sc.pl.umap(adata, color=["scissor_survival", "cell_type"], size=1)

# %% [markdown]
# Scissor+ cells are associated with late stage or with having the corresponding mutation. 

# %%
adata.obs["scissor_tumor_stage"] = [
    {"scissor+": "late_stage", "scissor-": "early_stage"}.get(x, np.nan)
    for x in adata.obs["scissor_tumor_stage"]
]

# %%
adata.obs["scissor_kras_mutation"] = [
    {"scissor+": "KRAS mut", "scissor-": "no KRAS mut"}.get(x, None)
    for x in adata.obs["scissor_kras_mutation"]
]
adata.obs["scissor_braf_mutation"] = [
    {"scissor+": "BRAF mut", "scissor-": "no BRAF mut"}.get(x, None)
    for x in adata.obs["scissor_braf_mutation"]
]
adata.obs["scissor_egfr_mutation"] = [
    {"scissor+": "EGFR mut", "scissor-": "no EGFR mut"}.get(x, None)
    for x in adata.obs["scissor_egfr_mutation"]
]

# %%
adata_primary = adata[adata.obs["origin"] == "tumor_primary", :].copy()

# %%
sc.pl.umap(adata_primary, color="scissor_tumor_stage", size=2)
for var in ["scissor_kras_mutation", "scissor_braf_mutation", "scissor_egfr_mutation"]:
    sc.pl.umap(adata_primary, color=var, size=2, palette=["#ff7f0e", "#1f77b4"])

# %%
sc.pl.umap(adata, color="VEGFA", size=1)

# %%
sc.pl.umap(adata, color=["PDCD1", "LAG3", "HAVCR2", "CXCL13", "CTLA4"], size=1), 

# %%
sc.pl.umap(adata, color=["dataset", "condition", "origin"], size=1)

# %%
adata_epi.obs["scissor"] = adata.obs["scissor_survival"]
adata_epi.obs["cell_type"] = adata.obs["cell_type"]
adata_tumor.obs["scissor"] = adata.obs["scissor_survival"]

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
pd.__version__


# %%
def scissor_by_group(adata, *, groupby="cell_type", scissor_col, adatas_for_gini=None):
    """Aggregate scissor scores first by patient, then by a grouping variable

    Parameters
    ----------
    adatas_for_gini
        Set this to a dict {groupby: adata} with an AnnData object
        for each group with a `leiden` clustering. This will
        be used to obtain gini index.
    """
    obs = adata_primary.obs.copy()
    # convert to str that nans are counted
    obs[scissor_col] = obs[scissor_col].astype(str)
    df_grouped = (
        (
            obs.groupby(["cell_type", "patient"], observed=True)[scissor_col]
            .value_counts(normalize=True)
            .reset_index(name="frac")
        )
        .pivot_table(
            values="frac",
            columns=scissor_col,
            index="cell_type",
            fill_value=0,
            aggfunc=np.mean,
        )
        .reset_index()
    )

    return df_grouped


# %%
def plot_scissor_df(df, *, title="scissor"):
    """Plot the result of scissor_by_group as a bar chart"""
    up, down = [x for x in df.columns[1:] if x != "nan"]
    order = df.sort_values(up)["cell_type"].values.tolist()
    return alt.vconcat(
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X("cell_type", sort=order, axis=None),
            y=up,
            # color=alt.Color(
            #     "gini_better", scale=alt.Scale(scheme="magma", reverse=False)
            # ),
        )
        .properties(height=100, title=title),
        alt.Chart(df.assign(**{down: lambda x: -x[down]}))
        .mark_bar()
        .encode(
            x=alt.X("cell_type", sort=order),
            y=down,
            # color=alt.Color(
            #     "gini_worse", scale=alt.Scale(scheme="magma", reverse=False)
            # ),
        )
        .properties(height=100),
        spacing=0,
    )


# %%
def subcluster_adata(tmp_adata):
    sc.pp.neighbors(tmp_adata, use_rep="X_scANVI")
    sc.tl.leiden(tmp_adata, resolution=0.5)
    return adata


# %%
adatas = []
cell_types = adata_primary.obs["cell_type"].unique()
for cell_type in cell_types:
    adatas.append(adata_primary[adata_primary.obs["cell_type"] == cell_type, :].copy())

# %%
adatas = process_map(subcluster_adata, adatas)


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
scissor_dfs = {
    k: scissor_by_group(adata, scissor_col=k)
    for k in [
        "scissor_survival",
        "scissor_tumor_stage",
        "scissor_kras_mutation",
        "scissor_braf_mutation",
        "scissor_egfr_mutation",
    ]
}

# %%
for col, df in scissor_dfs.items():
    plot_scissor_df(df, title=col).display()

# %% [markdown] tags=[]
# ---
#
# # Signature enrichment

# %%
sc.settings.set_figure_params(figsize=(4, 4))

# %%
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


# %% [markdown]
# ---

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
dfs = {}
for signature, adatas in all_adatas.items():
    tmp_dfs = []
    for ct, ad in zip(cell_types, adatas):
        try:
            tmp_dfs.append(paired_test(ad).assign(cell_type=ct))
        except ValueError:
            pass
    dfs[signature] = pd.concat(tmp_dfs)

# %%
tumor_normal_df = pd.concat(dfs).sort_values("pvalue")

# %%
_, tumor_normal_df["fdr"] = statsmodels.stats.multitest.fdrcorrection(
    tumor_normal_df["pvalue"].values,
    alpha=0.2,
)

# %%
pd.set_option("display.max_rows", None)

# %%
tumor_normal_df.shape

# %% [markdown]
# ---

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

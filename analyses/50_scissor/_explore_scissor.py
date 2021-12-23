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
scissor_res_files = {
    id: Path(
        "../../data/20_integrate_scrnaseq_data/scissor/scissor_by_sample/artifacts/"
    ).glob(f"*_scissor_{id}.tsv")
    for id in [
        "survival",
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

# %% [markdown]
# ---
#
# # Neutrophils

# %%
adata_neutro = adata[adata.obs["cell_type"] == "Granulocytes", :].copy()

# %%
sc.settings.set_figure_params(figsize=(5, 5))

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
    vmax=1
)

# %%
sc.pl.dotplot(
    neutro_ukimv[neutro_ukimv.obs["origin"] == "normal_adjacent", :],
    groupby="patient",
    var_names="VEGFA",
    title="normal_adjacent",
    vmin=0, 
    vmax=1
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
    vmax=1
)

# %%
adata_neutro.obs.groupby(
    ["dataset", "patient", "origin"], observed=True
).size().reset_index(name="n")

# %%
sc.pl.umap(adata_neutro, color=["origin", "dataset"])

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
    ad_ctf,
    var_names=ad_ctf.var_names,
    groupby="patient",
    dendrogram=True,
    swap_axes=True,
    cmap="bwr",
    vmin=-0.25,
    vmax=0.25
    # vmin=0,
    # vmax=1,
    # standard_scale="var",
)

# %%
sc.pp.neighbors(ad_ctf, use_rep="X")

# %%
sc.tl.umap(ad_ctf)

# %%
sc.pl.umap(
    ad_ctf, color=ad_ctf.var_names, cmap="bwr", vmin=-0.5, vmax=0.5, add_outline=True
)

# %%
sc.tl.leiden(ad_ctf)

# %%
sc.pl.umap(ad_ctf, color=["leiden", "dataset", "condition"])

# %%
sc.tl.rank_genes_groups(ad_ctf, groupby="leiden", method="wilcoxon")

# %%
sc.pl.rank_genes_groups_matrixplot(
    ad_ctf, values_to_plot="log10_pvals_adj", dendrogram=False, n_genes=4, vmax=5
)

# %%
sc.pl.matrixplot(
    ad_ctf,
    var_names=ad_ctf.var_names,
    groupby="leiden",
    cmap="bwr",
    vmax=0.25,
    vmin=-0.25,
    swap_axes=True,
)

# %%
# ad_ctf.obs.loc[:, ["patient_group", "patient"]].to_csv("patient_stratification.csv")

# %%
ah.annotate_cell_types(
    ad_ctf,
    {
        "T/B cell infiltration": [0, 5],
        "neutral": [1],
        "immune-excluded": [2, 3],
        "myeloid infiltration": [4],
        "Club": [6],
    },
    key_added="patient_group",
)

# %%
sc.pl.matrixplot(
    ad_ctf,
    var_names=ad_ctf.var_names,
    groupby="patient_group",
    cmap="bwr",
    vmax=0.25,
    vmin=-0.25,
    swap_axes=False,
)

# %%
patient_map = {p: g for p, g in zip(ad_ctf.obs["patient"], ad_ctf.obs["patient_group"])}

# %%
adata.obs["patient_group"] = [patient_map.get(p, None) for p in adata.obs["patient"]]

# %%
sc.pl.umap(adata, color=["patient_group", "cell_type"])

# %%
adata_epi.obs["patient_group"] = adata.obs["patient_group"]

# %%
with plt.rc_context({"figure.figsize": (7, 7)}):
    sc.pl.umap(adata_epi, color=["patient_group"], groups=["Club"], size=2)
    sc.pl.umap(adata_epi, color=["origin", "cell_type"], size=2)

# %% [markdown]
# ---
#
# # Signature enrichment

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
def run_dorothea(adata):
    tmp_adata = adata.copy()
    dorothea.run(
        tmp_adata,  # Data to use
        regulons,  # Dorothea network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # TF with less than 5 targets will be ignored
    )
    return dorothea.extract(tmp_adata)


# %%
def run_cytosig(adata):
    tmp_adata = adata.copy()
    progeny.run(
        tmp_adata,  # Data to use
        cytosig_signature,  # PROGENy network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # Pathways with less than 5 targets will be ignored
        obsm_key="cytosig",
    )
    return progeny.extract(tmp_adata, obsm_key="cytosig")


# %%
adatas_by_cell_type = []
cell_types = adata.obs["cell_type"].unique()
for cell_type in cell_types:
    adatas_by_cell_type.append(adata[adata.obs["cell_type"] == cell_type, :].copy())

# %%
adatas_dorothea = list(
    tqdm(map(run_dorothea, adatas_by_cell_type), total=len(adatas_by_cell_type))
)

# %% tags=[]
adatas_progeny = list(
    tqdm(map(run_progeny, adatas_by_cell_type), total=len(adatas_by_cell_type))
)

# %%
adatas_cytosig = list(
    tqdm(map(run_cytosig, adatas_by_cell_type), total=len(adatas_by_cell_type))
)

# %%
all_adatas = {
    "progeny": adatas_progeny,
    "cytosig": adatas_cytosig,
    "dorothea": adatas_dorothea,
}


# %%
# def test_lm(pseudobulk, groupby, covariate_formula, groups="all", reference="rest"):
#     """tmp_adata is a pseudobulk anndata object"""
#     var_names = pseudobulk.var_names
#     all_groups = set(pseudobulk.obs[groupby])
#     if groups == "all":
#         groups = all_groups

#     df = pseudobulk.obs.join(
#         pd.DataFrame(pseudobulk.X, index=pseudobulk.obs_names, columns=var_names)
#     )

#     pvalues = []
#     for col in var_names:
#         for group in groups:
#             group = set([group])
#             if reference == "rest":
#                 tmp_ref = all_groups - group
#             else:
#                 tmp_ref = set(reference)

#             tmp_df = df.loc[lambda x: x[groupby].isin(group | tmp_ref), :].copy()
#             tmp_df[groupby] = [
#                 "group" if x in group else "reference" for x in tmp_df[groupby]
#             ]

#             mod = smf.ols(formula=f"{col} ~ {groupby} {covariate_formula}", data=tmp_df)
#             res = mod.fit()

#             return res

#     # return pd.DataFrame(
#     #     index=var_names, data={"mean1": mean1, "mean2": mean2, "pvalue": pvalues}
#     # )

# %%
def test_lm(pseudobulk, formula):
    """
    Use a linear model to find differences between groups

    In this case we use sum-to-zero or deviation coding to find
    deviations from the mean of means

    tmp_adata is a pseudobulk anndata object"""
    var_names = pseudobulk.var_names

    df = pseudobulk.obs.join(
        pd.DataFrame(pseudobulk.X, index=pseudobulk.obs_names, columns=var_names)
    )

    pseudobulk.obs["patient_group"] = pd.Categorical(pseudobulk.obs["patient_group"])
    all_groups = pseudobulk.obs["patient_group"].unique()

    def test_all_params(res, all_groups):
        # only using the categories gets rid of NAN
        keys = [f"C(patient_group, Sum)[S.{g}]" for g in all_groups.categories.values]
        # print(keys)
        # print(res.params)
        coefs = res.params[keys[:-1]].to_dict()
        pvals = res.pvalues[keys[:-1]].to_dict()
        # test the level that was omitted for redundancy
        coefs[keys[-1]] = -sum(coefs.values())
        pvals[keys[-1]] = float(
            res.f_test(" + ".join([f"{k}" for k in keys[:-1]]) + " = 0").pvalue
        )
        return coefs, pvals

    results = []
    lms = []
    for col in var_names:
        group_results = []
        # there must be a better way to get all pvalues instead of re-training the LM for each group!
        # I can get n-1 pvalues from the model, there must be a way to get the last one!
        mod = smf.ols(formula=formula.format(col=col), data=df)
        res = mod.fit()
        coefs, pvals = test_all_params(res, all_groups)
        res_df = (
            pd.DataFrame.from_dict(coefs, orient="index", columns=["coef"])
            .join(pd.DataFrame.from_dict(pvals, orient="index", columns=["pvalue"]))
            .assign(
                variable=col,
                group=lambda x: [
                    re.search("\[S\.(.*)\]", k).groups()[0] for k in x.index
                ],
            )
        )
        results.append(res_df)
        lms.append(res)

    return lms, pd.concat(results)


# %% tags=[] jupyter={"outputs_hidden": true}
res_per_cell_type = {}
for signature, adatas in all_adatas.items():
    tmp_res = []
    for ct, tmp_adata in zip(tqdm(cell_types), adatas):
        tmp_bdata = hb.tl.pseudobulk(
            tmp_adata, groupby=["dataset", "patient", "patient_group"], aggr_fun=np.mean
        )
        if tmp_bdata.obs["patient_group"].nunique() < 3:
            continue
        _, res_df = test_lm(
            tmp_bdata, "Q('{col}') ~ 0 + C(patient_group, Sum) + dataset"
        )
        res_df = res_df.reset_index(drop=True).assign(cell_type=ct)
        tmp_res.append(res_df)

    res_per_cell_type[signature] = pd.concat(tmp_res).dropna(how="any")

# %%
for tmp_df in res_per_cell_type.values():
    _, tmp_df["fdr"] = statsmodels.stats.multitest.fdrcorrection(
        tmp_df["pvalue"].values
    )

# %%
res_per_cell_type["progeny"].sort_values("pvalue").groupby("group").apply(
    lambda x: x.head(10)
)

# %%
d_progeny = {ct: ad for ct, ad in zip(cell_types, adatas_progeny)}


# %%

# %%
def mk_matrixplot(cell_type, features):
    tmp_bulk = hb.tl.pseudobulk(
        d_progeny[cell_type], groupby=["patient", "patient_group"], aggr_fun=np.mean
    )
    sc.pl.matrixplot(
        tmp_bulk, var_names=features, groupby="patient_group", title=cell_type
    )


# %%
mk_matrixplot("T cell CD4", ["PI3K", "VEGF", "TNFa"])

# %%
mk_matrixplot("Macrophage", ["p53", "PI3K", "Hypoxia", "TGFb"])

# %%
res_per_cell_type["dorothea"].sort_values("pvalue").groupby("group").apply(
    lambda x: x.head()
)

# %% tags=[]
res_per_cell_type["cytosig"].sort_values("pvalue").groupby("group").apply(
    lambda x: x.head()
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

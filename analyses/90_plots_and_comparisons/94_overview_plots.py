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
from scanpy_helpers.annotation import AnnotationHelper
from nxfvars import nxfvars
import altair as alt
from toolz.functoolz import pipe
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy_helpers as sh
import statsmodels.formula.api as smf

# %%
alt.data_transformers.disable_max_rows()

# %%
sh.colors.plot_all_palettes()

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
main_adata = nxfvars.get(
    "main_adata",
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
)

# %%
adata = sc.read_h5ad(main_adata)

# %%
with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": 300}):
    sc.pl.umap(
        adata,
        color="cell_type_coarse",
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=3,
    )

# %%
adata.shape

# %%
adata.obs["sample"].nunique()

# %%
adata.obs["platform"].nunique()

# %%
adata.obs["cell_type"].nunique()

# %%
df = (
    adata.obs.loc[lambda x: x["origin"].str.contains("tumor")]
    .groupby(["dataset"])
    .apply(lambda x: x["patient"].nunique())
    .reset_index(name="# lung cancer patients")
    .loc[lambda x: x["# lung cancer patients"] > 0]
    .sort_values("# lung cancer patients")
)
ch = (
    alt.Chart(df)
    .mark_bar()
    .encode(
        x="# lung cancer patients",
        y=alt.Y("dataset", sort="x"),
        color=alt.Color(
            "# lung cancer patients", scale=alt.Scale(scheme="viridis"), legend=None
        ),
    )
)
ch + (ch).mark_text(align="left", dx=3).encode(
    text="# lung cancer patients", color=alt.value("black")
)


# %%
def process_subset(mask):
    print(np.sum(mask))
    adata_sub = adata[mask, :].copy()
    sc.pp.neighbors(adata_sub, use_rep="X_scANVI")
    sc.tl.umap(adata_sub)
    return adata_sub


# %%
adatas = {
    label: process_subset(ct)
    for label, ct in {
        "epithelial": (adata.obs["cell_type_coarse"] == "Epithelial cell")
        & (adata.obs["cell_type_major"] != "other"),
        "tumor": (adata.obs["cell_type"] == "Tumor cells"),
        "immune": adata.obs["cell_type_coarse"].isin(
            [
                "T cell",
                "NK cell",
                "Neutrophils",
                "Plasma cell",
                "Mast cell",
                "B cell",
                "pDC",
                "cDC",
                "Macrophage/Monocyte",
            ]
        ),
        "structural": adata.obs["cell_type_coarse"].isin(
            ["Endothelial cell", "Stromal"]
        ),
    }.items()
}

# %%
with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": 300}):
    sc.pl.umap(
        adata,
        color="cell_type_coarse",
        legend_loc=None,
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=3,
        groups=["Epithelial cell"],
    )

# %%
with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": 300}):
    tmp_adata = adatas["epithelial"].copy()
    tmp_adata.obs["cell_type_tumor"] = tmp_adata.obs["cell_type_tumor"].str.replace(
        "Tumor cells ", ""
    )
    fig = sc.pl.umap(
        tmp_adata,
        color="cell_type_tumor",
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=8,
        return_fig=True,
    )
    fig.savefig("/home/sturm/Downloads/umap_epi.pdf", dpi=600)
    plt.show()

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(
        adatas["tumor"],
        color="cell_type_tumor",
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        add_outline=True,
        size=12,
    )

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(
        adatas["immune"],
        color="cell_type_coarse",
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        add_outline=True,
        size=12,
    )

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(
        adatas["immune"],
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        add_outline=True,
        size=12,
    )

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(
        adatas["structural"],
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        add_outline=True,
        size=12,
    )

# %% [markdown]
# ### By cell type

# %%
count_by_cell_type = adata.obs.groupby(["dataset", "cell_type_coarse"]).agg(
    n_cells=pd.NamedAgg("cell_type_coarse", "count")
)
count_by_cell_type["frac_cells"] = adata.obs.groupby("dataset")[
    "cell_type_coarse"
].value_counts(normalize=True)
count_by_cell_type.reset_index(inplace=True)
count_by_cell_type.to_csv(f"{artifact_dir}/cell_type_coarse_per_dataset.tsv", sep="\t")

# %%
c1 = (
    alt.Chart(count_by_cell_type)
    .mark_bar()
    .encode(x="n_cells", color="cell_type_coarse", y="dataset")
)
s1 = (
    alt.Chart(count_by_cell_type)
    .mark_bar()
    .encode(x="n_cells", color="cell_type_coarse")
)
c1


# %% [markdown]
# ### By sample

# %%
data = (
    adata.obs.groupby(["dataset", "sample", "origin"], observed=True)
    .size()
    .reset_index(name="n_cells")
)
c2 = (
    alt.Chart(data)
    .mark_bar()
    .encode(
        x=alt.X("count(sample)", title="n_samples"),
        color="origin",
        y=alt.Y("dataset", axis=None),
    )
)
s2 = (
    alt.Chart(data)
    .mark_bar()
    .encode(x=alt.X("count(sample)", title="n_samples"), color="origin")
)
c2

# %% [markdown]
# ### By patients

# %%
data = (
    adata.obs.groupby(["dataset", "patient", "condition"], observed=True)
    .size()
    .reset_index(name="n_cells")
)
c3 = (
    alt.Chart(data)
    .mark_bar()
    .encode(
        x=alt.X("count(patient)", title="n_patients"),
        color="condition",
        y=alt.Y("dataset", axis=None),
    )
)
s3 = (
    alt.Chart(data)
    .mark_bar()
    .encode(x=alt.X("count(patient)", title="n_patients"), color="condition")
)
c3


# %%
(
    s1.properties(width=200) | s2.properties(width=200) | s3.properties(width=200)
).resolve_scale(color="independent")

# %%
(
    c1.properties(width=200) | c2.properties(width=200) | c3.properties(width=200)
).resolve_scale(color="independent")

# %% [markdown]
# ### cell-type composition (relative) by dataset

# %%
c1.encode(x=alt.X("n_cells", title="n_cells", stack="normalize"))

# %% [markdown]
# ### Comparison of cell-type composition between sample origin

# %%
adata.obs.groupby(["origin", "cell_type_coarse"]).size().reset_index(
    name="n_cells"
).loc[
    lambda x: x["origin"].isin(["normal", "normal_adjacent", "tumor_primary"]), :
].pipe(
    lambda x: alt.Chart(x)
    .mark_bar()
    .encode(x=alt.X("n_cells", stack="normalize"), y="origin", color="cell_type_coarse")
)

# %%
adata.obs.loc[
    lambda x: (x["origin"] == "tumor_primary")
    & ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])
].groupby(["condition", "cell_type_coarse"]).size().reset_index(name="n_cells").loc[
    lambda x: x["condition"].isin(["LUAD", "LSCC"]), :
].pipe(
    lambda x: alt.Chart(x)
    .mark_bar()
    .encode(
        x=alt.X("n_cells", stack="normalize"),
        y="condition",
        color=alt.Color(
            "cell_type_coarse", scale=sh.colors.altair_scale("cell_type_coarse")
        ),
    )
)

# %%
adata.obs.loc[
    lambda x: (x["origin"] == "tumor_primary")
    & ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])
].groupby(["tumor_stage", "cell_type_coarse"]).size().reset_index(name="n_cells").pipe(
    lambda x: alt.Chart(x)
    .mark_bar()
    .encode(
        x=alt.X("n_cells", stack="normalize"),
        y="tumor_stage",
        color=alt.Color(
            "cell_type_coarse", scale=sh.colors.altair_scale("cell_type_coarse")
        ),
    )
)

# %% [markdown]
# ### Compare cell-type fractions between LUAD und LSCC

# %%
import sccoda.util.cell_composition_data as scc_dat
import sccoda.util.comp_ana as scc_ana
import sccoda.util.data_visualization as scc_viz

# %%
frac_by_condition = (
    adata.obs.loc[
        lambda x: (x["origin"] == "tumor_primary")
        & ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])
        & x["condition"].isin(["LUAD", "LSCC"]) # & 
        # (x["cell_type_coarse"] != "Epithelial cell")
    ]
    .groupby(["dataset", "condition", "tumor_stage", "patient"])
    .apply(lambda x: x.value_counts("cell_type_major", normalize=False))
    .reset_index(name="n_cells")
    .assign(condition=lambda x: x["condition"].astype(str))
)

# %%
frac_by_condition

# %%
frac_pivot = frac_by_condition.pivot(
    index=["patient", "dataset", "condition", "tumor_stage"], columns="cell_type_major", values="n_cells"
).reset_index()

# %%
data_all = scc_dat.from_pandas(frac_pivot, covariate_columns=["patient", "dataset", "condition", "tumor_stage"])

# %%
sccoda_mod = scc_ana.CompositionalAnalysis(data_all, formula="condition + dataset", reference_cell_type="Tumor cells")

# %%
sscoda_res = sccoda_mod.sample_hmc(num_results = 10000)

# %%
scc_viz.stacked_barplot(data_all, feature_name="condition")

# %%
scc_viz.boxplots(data_all, feature_name="condition", figsize=(12, 5))

# %%
sscoda_res.summary()

# %%
sscoda_res.credible_effects(est_fdr=0.1).reset_index().head(20)

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 8))
sns.boxplot(
    x="frac",
    y="cell_type_coarse",
    hue="condition",
    data=frac_by_condition,
    palette=sh.colors.COLORS.condition,
    fliersize=0,
    ax=ax,
)
sns.stripplot(
    x="frac",
    y="cell_type_coarse",
    hue="condition",
    data=frac_by_condition,
    ax=ax,
    dodge=True,
    color="black",
    size=3,
    alpha=0.5,
)

# %%
results = []
for ct in adata.obs["cell_type_coarse"].unique():
    mod = smf.ols(
        f"Q('{ct}') ~ C(condition) + dataset",
        data=frac_pivot,
    )
    res = mod.fit()
    results.append(
        {
            "cell_type": ct,
            "pvalue": res.pvalues["C(condition)[T.LUAD]"],
            "coef": res.params["C(condition)[T.LUAD]"],
        }
    )

# %%
pd.DataFrame(results).sort_values("pvalue")

# %% [markdown]
# ### Compare cell-type fractions between early and late

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 8))
sns.boxplot(
    x="frac",
    y="cell_type_coarse",
    hue="tumor_stage",
    data=frac_by_condition,
    palette=sh.colors.COLORS.tumor_stage,
    fliersize=0,
    ax=ax,
)
sns.stripplot(
    x="frac",
    y="cell_type_coarse",
    hue="tumor_stage",
    data=frac_by_condition,
    ax=ax,
    dodge=True,
    color="black",
    size=3,
    alpha=0.5,
)

# %%
frac_pivot = frac_by_condition.pivot(
    index=["patient", "dataset", "tumor_stage"], columns="cell_type_coarse", values="frac"
).reset_index()

# %%
results = []
for ct in adata.obs["cell_type_coarse"].unique():
    mod = smf.ols(
        f"Q('{ct}') ~ C(tumor_stage) + dataset",
        data=frac_pivot,
    )
    res = mod.fit()
    results.append(
        {
            "cell_type": ct,
            "pvalue": res.pvalues["C(tumor_stage)[T.late]"],
            "coef": res.params["C(tumor_stage)[T.late]"],
        }
    )

# %%
pd.DataFrame(results).sort_values("pvalue")

# %% [markdown]
# ---

# %% [markdown]
# ## Plots by sample

# %%
sample_df = adata.obs.loc[
    :, ["sample", "tissue", "origin", "condition", "dataset"]
].drop_duplicates()
assert sample_df.shape[0] == adata.obs["sample"].nunique()

# %%
(
    sample_df.groupby("dataset")
    .size()
    .reset_index(name="n")
    .pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="dataset"))
)

# %%
(
    sample_df.groupby("tissue")
    .size()
    .reset_index(name="n")
    .pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="tissue"))
)

# %%
(
    sample_df.groupby("condition")
    .size()
    .reset_index(name="n")
    .pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="condition"))
)

# %%
origins = {ds: set() for ds in adata.obs["dataset"].unique()}

# %%
for ds, origin in zip(adata.obs["dataset"], adata.obs["origin"]):
    origins[ds] |= {origin}

# %%
for ds, origin in origins.items():
    if "normal" in origin and any([s.startswith("tumor") for s in origin]):
        print(ds, origin)

# %%
(
    sample_df.groupby("origin")
    .size()
    .reset_index(name="n")
    .pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="origin"))
)


# %% [markdown]
# ## Plots by patient

# %%
def _group_origins(origins):
    origins = set(origins)
    if len(origins) == 1:
        return next(iter(origins))
    elif "tumor_primary" in origins and "normal_adjacent" in origins:
        return "matched tumor/normal"
    else:
        # print(origins)
        return "multiple"


patient_df = adata.obs.groupby(["patient", "dataset", "condition"], observed=True).agg(
    {"origin": _group_origins}
)

assert patient_df.shape[0] == adata.obs["patient"].nunique()

# %%
patient_df

# %%
(
    patient_df.groupby("origin")
    .size()
    .reset_index(name="n")
    .pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="origin"))
)

# %%
(
    patient_df.groupby("dataset")
    .size()
    .reset_index(name="n")
    .pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="dataset"))
)

# %%
(
    patient_df.groupby("condition")
    .size()
    .reset_index(name="n")
    .pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="condition"))
)

# %% [markdown]
# ---

# %%
exit

# %%
adata.obs["cell_type_plot"] = [
    "tumor cell" if x == "Tumor cells" else "other cell" for x in adata.obs["cell_type"]
]

# %%
composition_df = (
    adata.obs.groupby(["patient", "cell_type"]).size().reset_index(name="n")
)

# %%
patient_order = (
    composition_df
    >> s.group_by("patient")
    >> s.mutate(frac=_.n / _.n.sum())
    >> s.filter(_.cell_type == "Tumor cells")
    >> s.ungroup()
    >> s.arrange("frac")
)["patient"].values.tolist()

# %%
bar_chart = (
    alt.Chart(composition_df)
    .mark_bar()
    .encode(
        x=alt.X(
            "patient", sort=patient_order, axis=alt.Axis(labels=False, ticks=False)
        ),
        color=alt.Color("cell_type", legend=alt.Legend(columns=2, symbolLimit=0)),
        y=alt.Y("sum(n)", stack="normalize"),
    )
    .properties(height=400, width=1200)
)

bar_chart

# %%
condition_annotation = (
    alt.Chart(patient_df)
    .mark_rect()
    .encode(
        x=alt.X(
            "patient:N", sort=patient_order, axis=alt.Axis(labels=False, ticks=False)
        ),
        color="condition",
    )
    .properties(width=1200)
)
dataset_annotation = (
    alt.Chart(patient_df)
    .mark_rect()
    .encode(
        x=alt.X(
            "patient:N", sort=patient_order, axis=alt.Axis(labels=False, ticks=False)
        ),
        color="dataset",
    )
    .properties(width=1200)
)

origin_annotation = (
    alt.Chart(patient_df)
    .mark_rect()
    .encode(
        x=alt.X(
            "patient:N", sort=patient_order, axis=alt.Axis(labels=False, ticks=False)
        ),
        color="origin",
    )
    .properties(width=1200)
)


(condition_annotation & dataset_annotation & origin_annotation & bar_chart)

# %% [markdown]
# ## tumor cells only

# %%
obs_tumor = adata_tumor.obs.loc[
    adata_tumor.obs.origin.str.startswith("tumor"), :
].reset_index()

# %%
obs_tumor.head()

# %%
tumor_by_cell_type = (
    obs_tumor.groupby(["patient", "cell_type"], observed=True)
    .size()
    .reset_index(name="n")
)

# %%
tumor_sort_patients = (
    tumor_by_cell_type.groupby("patient", observed=True)
    .apply(
        lambda x: x.assign(
            n_luad=x.loc[x["cell_type"].str.contains("LUAD"), "n"].sum(), n=x.n.sum()
        )
    )
    .drop("cell_type", axis=1)
    .drop_duplicates()
    .assign(luad_frac=lambda x: x["n_luad"] / x["n"])
    .sort_values("luad_frac")["patient"]
    .values.tolist()
)

# %%
bar_chart = (
    alt.Chart(tumor_by_cell_type)
    .mark_bar()
    .encode(
        x=alt.X(
            "patient",
            sort=tumor_sort_patients,
            axis=alt.Axis(labels=False, ticks=False),
        ),
        color=alt.Color("cell_type", legend=alt.Legend(columns=2, symbolLimit=0)),
        y=alt.Y("sum(n)", stack="normalize"),
    )
    .properties(height=400, width=1200)
)

# %%
condition_annotation = (
    alt.Chart(
        obs_tumor.set_index("patient")
        .loc[tumor_sort_patients, ["condition"]]
        .reset_index()
        .drop_duplicates()
    )
    .mark_rect()
    .encode(
        x=alt.X(
            "patient:N",
            sort=tumor_sort_patients,
            axis=alt.Axis(labels=False, ticks=False),
        ),
        color="condition",
    )
    .properties(width=1200)
)

bar_chart & condition_annotation

# %%

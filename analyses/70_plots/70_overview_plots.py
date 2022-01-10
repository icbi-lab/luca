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
import scanpy as sc
from scanpy_helpers.annotation import AnnotationHelper
from nxfvars import nxfvars
import altair as alt
from toolz.functoolz import pipe
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %%
alt.data_transformers.disable_max_rows()

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")

# %%
main_adata = nxfvars.get(
    "main_adata",
    "../../data/20_build_atlas//annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)

# %%
adata = sc.read_h5ad(main_adata)

# %%
adata_epi = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/epithelial_cells_annotated.h5ad"
)

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(
        adata,
        color="cell_type_coarse",
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        add_outline=True,
        size=2,
    )

# %%
adata.obs["cell_type"].value_counts()

# %%
with plt.rc_context({"figure.figsize": (8, 8)}):
    sc.pl.umap(
        adata_epi,
        color="cell_type",
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        frameon=False,
        add_outline=True,
        size=4,
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
s1 = alt.Chart(data).mark_bar().encode(x="n_cells", color="cell_type_coarse")
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

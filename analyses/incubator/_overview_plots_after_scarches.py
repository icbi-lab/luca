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
from tqdm.contrib.concurrent import process_map
import itertools
import progeny
import dorothea
from threadpoolctl import threadpool_limits

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
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
)
epithelial_adata = nxfvars.get(
    "main_adata",
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/epithelial_cells_annotated.h5ad",
)
threadpool_limits(nxfvars.get("cpus", 16))

# %%
adata = sc.read_h5ad(main_adata)

# %%
adata_epi = sc.read_h5ad(epithelial_adata)

# %% [markdown]
# # Stats and metadata

# %%
## Export patient table
patient_metadata = adata.obs.loc[
    :,
    [
        # "study",
        "dataset",
        "patient",
        "uicc_stage",
        "tumor_stage",
        "sex",
        "ever_smoker",
        "driver_genes",
        "condition",
        "age",
        "platform",
        "platform_fine",
    ] + [x for x in adata.obs.columns if "mutation" in x],
].drop_duplicates().sort_values(
    [
        # "study",
        "dataset",
        "patient",
    ]
).reset_index(
    drop=True
)

patient_metadata.to_csv(
    f"{artifact_dir}/patient_table.csv"
)

# %%
pd.set_option("display.max_rows", 500)

# %%
patient_metadata

# %% [markdown]
# # UMAP plots
#
# UMAP overview plots of the "core" atlas (without projected datasets)

# %%
with plt.rc_context({"figure.figsize": (5, 5), "figure.dpi": 300}):
    fig = sc.pl.umap(
        adata,
        color="cell_type_coarse",
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontoutline=2,
        frameon=False,
        # add_outline=True,
        size=1,
        return_fig=True,
        title="",
    )
    fig.savefig(f"{artifact_dir}/umap_core_atlas.pdf")

# %% [markdown]
# # Stats

# %%
adata.obs["sample"].nunique()

# %%
adata.obs["platform"].nunique()

# %%
adata.obs["cell_type"].nunique()

# %%
adata.obs["patient"].nunique()

# %%
adata.obs["study"].nunique()

# %%
adata.obs.columns

# %%
adata.obs.loc[lambda x: x["origin"].str.contains("tumor")]["patient"].nunique()

# %%
df = (
    adata.obs.loc[lambda x: x["origin"].str.contains("tumor")]
    .groupby(["study"])
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
        y=alt.Y("study", sort="x"),
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
        "tumor": (adata.obs["cell_type_major"] == "Tumor cells")
        & (adata.obs["cell_type_tumor"] != "Club"),
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

# %% [markdown]
# ### Tumor stages

# %%
adata.obs.loc[:, ["dataset", "patient", "tumor_stage"]].drop_duplicates().groupby(
    ["dataset", "tumor_stage"]
).size().unstack()

# %% [markdown]
# ---

# %% [markdown]
# ### Progeny and dorothea for tumor cell clusters

# %%
model = progeny.load_model(
    organism="Human",  # If working with mouse, set to Mouse
    top=1000,  # For sc we recommend ~1k target genes since there are dropouts
)

# %%
progeny.run(
    adatas["tumor"],  # Data to use
    model,  # PROGENy network
    center=True,  # Center gene expression by mean per cell
    num_perm=100,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # Pathways with less than 5 targets will be ignored
)

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
    adatas["tumor"],  # Data to use
    regulons,  # Dorothea network
    center=True,  # Center gene expression by mean per cell
    num_perm=100,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
    use_raw=True,  # Use raw adata, where we have the lognorm gene expression
    min_size=5,  # TF with less than 5 targets will be ignored
)

# %%
ad_tumor_progeny = progeny.extract(adatas["tumor"])

# %%
sc.pl.matrixplot(
    ad_tumor_progeny,
    var_names=ad_tumor_progeny.var_names,
    groupby="cell_type_tumor",
    cmap="bwr",
    vmin=-1,
    vmax=1,
)

# %%
ad_tumor_dorothea = dorothea.extract(adatas["tumor"])

# %%
sc.pl.matrixplot(
    ad_tumor_dorothea,
    var_names=ad_tumor_dorothea.var_names,
    groupby="cell_type_tumor",
    cmap="bwr",
    vmin=-3,
    vmax=3,
)

# %% [markdown]
# # UMAP plots

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
    lambda x: x["condition"].isin(["LUAD", "LUSC"]), :
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

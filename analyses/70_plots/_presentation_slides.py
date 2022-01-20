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
# %config InlineBackend.figure_formats = ['png']
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
import threadpoolctl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st
import matplotlib
import scanpy_helpers as sh
import altair as alt

# %%
threadpoolctl.threadpool_limits(20)

# %%
matplotlib.rcParams.update({"font.size": 16})
matplotlib.rcParams["figure.dpi"] = 72

# %%
adata = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad"
)

# %% [markdown]
# ## Cell types by platform

# %%
df = (
    adata.obs.loc[lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])]
    .groupby(["platform", "dataset"], observed=True)["cell_type_coarse"]
    .value_counts(normalize=True)
    .reset_index()
    .rename(columns={"cell_type_coarse": "fraction", "level_2": "cell_type_coarse"})
)

# %%
df

# %%
fig, ax = plt.subplots(1, 1, figsize=(7, 12))
sns.barplot(
    data=df,
    y="cell_type_coarse",
    x="fraction",
    hue="platform",
    ax=ax,
    palette=sh.colors.COLORS.platform,
    ci=None,
    orient="h",
)
sns.stripplot(
    data=df,
    y="cell_type_coarse",
    x="fraction",
    hue="platform",
    ax=ax,
    dodge=True,
    edgecolor="black",
    orient="h",
    linewidth=1,
    palette=sh.colors.COLORS.platform,
)
ax.grid(axis="x")
ax.legend(
    *zip(
        *[
            (h, l)
            for h, l in zip(*ax.get_legend_handles_labels())
            if "Bar" in str(type(h))
        ]
    ),
    loc="center left",
    bbox_to_anchor=(1, 0.5)
)
plt.show()

# %% [markdown]
# ## Pie charts atlas composition

# %%
df = (
    (
        adata.obs.groupby(["dataset", "sample", "origin"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .assign(
            origin=lambda x: [
                {
                    "normal": "normal",
                    "normal_adjacent": "adj. normal",
                    "tumor_edge": "primary tumor",
                    "tumor_metastasis": "metastases",
                    "tumor_middle": "primary tumor",
                    "tumor_primary": "primary tumor",
                }.get(x, "other")
                for x in x["origin"]
            ]
        )
        .groupby("origin")
        .size()
        .reset_index(name="n_samples")
    )
    .set_index("origin")
    .loc[["normal", "adj. normal", "primary tumor", "metastases", "other"]]
)
fig, ax = plt.subplots()
ax.pie(
    df["n_samples"],
    labels=df.index,
    colors=[sh.colors.COLORS.origin[c] for c in df.index],
)
plt.show()

# %%
df = (
    adata.obs.groupby(["dataset", "patient", "condition"], observed=True)
    .size()
    .reset_index(name="n_cells")
    .assign(
        condition=lambda x: [
            x if x in ["COPD", "LUAD", "LSCC", "NSCLC", "healthy_control"] else "other"
            for x in x["condition"]
        ]
    )
    .groupby("condition")
    .size()
    .reset_index(name="n_patients")
)
fig, ax = plt.subplots()
ax.pie(
    df["n_patients"],
    labels=df["condition"],
    colors=[sh.colors.COLORS.condition[c] for c in df["condition"]],
)
plt.show()

# %%
df = (
    adata.obs.groupby(["dataset", "platform"], observed=True)
    .size()
    .reset_index(name="n_cells")
    .groupby("platform")
    .size()
    .reset_index(name="n_datasets")
).sort_values("n_datasets", ascending=False)


def fmt_label(pct):
    n = int(pct / 100 * sum(df["n_datasets"]))
    return n


fig, ax = plt.subplots()
ax.pie(
    df["n_datasets"],
    labels=df["platform"],
    autopct=fmt_label,
    colors=[sh.colors.COLORS.platform[c] for c in df["platform"]],
)
plt.show()

# %% [markdown]
# ## Pseudoreplication bias

# %%
neutro = adata[
    (adata.obs["cell_type_major"] == "Neutrophils")
    & (adata.obs["dataset"] == "UKIM-V"),
    :,
].copy()

# %%
neutro.obs["origin"] = [
    {"normal_adjacent": "adj. normal", "tumor_primary": "tumor"}[x]
    for x in neutro.obs["origin"]
]

# %%
neutro.obs["patient"] = [p.split("_")[1] for p in neutro.obs["patient"]]

# %%
df = (
    pd.DataFrame()
    .assign(
        VEGFA=neutro.raw[:, "VEGFA"].X.todense().A1,
        patient=neutro.obs["patient"].values.tolist(),
        origin=neutro.obs["origin"].values.tolist(),
    )
    .sort_values(["patient", "origin"])
)

# %% [markdown]
# ### using only P2

# %%
tmp_df = df.loc[df["patient"] == "P2"]

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
sns.boxplot(
    data=tmp_df,
    x="origin",
    y="VEGFA",
    ax=ax,
    width=0.5,
    # cut=0,
    # inner="box",
    color="lightgrey",
)
# sns.boxplot(data=df, x="origin", y="VEGFA", ax=ax)
sns.stripplot(
    data=tmp_df,
    x="origin",
    y="VEGFA",
    ax=ax,
    color=sns.palettes.color_palette("colorblind")[1],
    size=2,
    # palette=,
)
ax.set_title("VEGFA expression")
ax.set_ylabel("log norm counts")
plt.grid(axis="y")

# %%
st.mannwhitneyu(
    tmp_df.loc[lambda x: x["origin"] == "adj. normal", "VEGFA"],
    tmp_df.loc[lambda x: x["origin"] == "tumor", "VEGFA"],
)

# %% [markdown]
# ### all patients

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
sns.boxplot(
    data=df,
    x="origin",
    y="VEGFA",
    ax=ax,
    width=0.5,
    # cut=0,
    # inner="box",
    color="lightgrey",
)
# sns.boxplot(data=df, x="origin", y="VEGFA", ax=ax)
sns.stripplot(
    data=df, x="origin", y="VEGFA", ax=ax, size=2, hue="patient", palette="colorblind"
)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax.set_title("VEGFA expression")
ax.set_ylabel("log norm counts")
plt.grid(axis="y")

# %%
st.mannwhitneyu(
    df.loc[lambda x: x["origin"] == "adj. normal", "VEGFA"],
    df.loc[lambda x: x["origin"] == "tumor", "VEGFA"],
)

# %%
neutro_pb = sh.pseudobulk.pseudobulk(neutro, groupby=["patient", "origin"])

# %%
sc.pp.normalize_total(neutro_pb)
sc.pp.log1p(neutro_pb)

# %%
df2 = (
    pd.DataFrame()
    .assign(
        VEGFA=neutro_pb[:, "VEGFA"].X[:, 0],
        patient=neutro_pb.obs["patient"].values,
        origin=neutro_pb.obs["origin"].values,
    )
    .sort_values("patient")
)

# %%
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
sns.boxplot(data=df2, x="origin", y="VEGFA", ax=ax, width=0.5, color="lightgrey")
# sns.boxplot(data=df, x="origin", y="VEGFA", ax=ax)
sns.stripplot(
    data=df2, x="origin", y="VEGFA", ax=ax, size=20, hue="patient", palette="colorblind"
)
sns.lineplot(
    data=df2,
    x="origin",
    y="VEGFA",
    ax=ax,
    hue="patient",
    palette="colorblind",
    legend=False,
)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
ax.set_title("VEGFA expression")
ax.set_ylabel("log norm counts")
plt.grid(axis="y")

# %%
neutro.obs.groupby(["patient", "origin"]).size().reset_index(name="n_cells")

# %%
with plt.rc_context({"figure.figsize": (3, 3)}):
    sns.barplot(
        # data=neutro_pb.obs.loc[lambda x: x["origin"] == "adj. normal"],
        data=neutro_pb.obs,
        x="patient",
        y="n_obs",
        ci=None,
        palette="colorblind",
        order=sorted(neutro_pb.obs["patient"].unique()),
    )

# %%
st.mannwhitneyu(
    df2.loc[lambda x: x["origin"] == "normal_adjacent", "VEGFA"],
    df2.loc[lambda x: x["origin"] == "tumor_primary", "VEGFA"],
)

# %%

# %%

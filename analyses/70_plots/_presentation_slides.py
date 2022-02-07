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
import numpy as np

# %%
threadpoolctl.threadpool_limits(20)

# %%
matplotlib.rcParams.update({"font.size": 16})
matplotlib.rcParams["figure.dpi"] = 300

# %%
adata = sc.read_h5ad(
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad"
)

# %% [markdown]
# ## TCGA

# %%
tcga_clin = pd.read_csv("../../tables/tcga/clinical_data_for_scissor.tsv", sep="\t")
tcga_clin["type"] = tcga_clin["type"].str.replace("LUSC", "LSCC")

# %%
df = tcga_clin.groupby("type").size().reset_index(name="n")
ch = (
    alt.Chart(df)
    .mark_bar(size=40)
    .encode(
        x=alt.X("type"),
        y=alt.Y("n", axis=alt.Axis(labels=False, ticks=False)),
        color=alt.Color("type", scale=sh.colors.altair_scale("condition"), legend=None),
    )
)
(
    ch + (ch).mark_text(dx=0, dy=15, size=15).encode(text="n", color=alt.value("black"))
).properties(width=100)

# %%
df = (
    tcga_clin.loc[:, ["kras_mutation", "braf_mutation", "egfr_mutation"]]
    .sum()
    .reset_index(name="n")
    .rename(columns={"index": "mutation", "n": "frac"})
)
df["mutation"] = df["mutation"].str.replace("_mutation", "").str.upper()
df["frac"] = df["frac"] / tcga_clin.shape[0]

# %%
ch = alt.Chart(df).mark_bar().encode(x=alt.X("frac", axis=alt.Axis(labels=False, ticks=False)), y="mutation")
(
    ch
    + (ch)
    .mark_text(align="left", dx=3)
    .encode(
        text=alt.Text("frac", format=",.2f", ),
        color=alt.value("black"),
    )
).configure_view(strokeOpacity=0).configure_axis(grid=False).properties(width=150)

# %%
tcga_clin["kras_mutation"].sum()

# %% [markdown]
# ## Cell types by platform

# %%
df = (
    adata.obs.loc[lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])][
        "cell_type_coarse"
    ]
    .value_counts(normalize=True)
    .reset_index(name="fraction")
    .rename(columns={"index": "cell_type"})
)
alt.Chart(df).mark_bar().encode(
    color=alt.Color(
        "cell_type", scale=sh.colors.altair_scale("cell_type_coarse"), legend=None
    ),
    x="fraction",
    y=alt.Y("cell_type", sort="-x"),
)

# %%

# %%
df = (
    (
        adata.obs.loc[
            lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])
        ]["platform"].value_counts(normalize=True)
    )
    .reset_index(name="fraction")
    .rename(columns={"index": "platform"})
    .assign(y="all cells")
)

c_vline = (
    alt.Chart(df.loc[lambda x: x["platform"] == "10x"])
    .mark_rule(color="black", strokeWidth=1, strokeDash=[3, 3])
    .encode(x="fraction")
)


c1 = (
    alt.Chart(df)
    .mark_bar()
    .encode(
        x=alt.X(
            "fraction",
            axis=None,
        ),
        order=alt.Order("fraction", sort="descending"),
        color=alt.Color(
            "platform",
            scale=sh.colors.altair_scale("platform"),
            sort=alt.EncodingSortField("fraction"),
        ),
        y=alt.Y("y", title=None),
    )
) + c_vline

df = (
    adata.obs.loc[lambda x: ~x["dataset"].isin(["Guo_Zhang_2018", "Maier_Merad_2020"])]
    .groupby("cell_type_coarse")["platform"]
    .value_counts(normalize=True)
    .reset_index(name="fraction")
    .rename(columns={"level_1": "platform"})
    .set_index("platform")
    .assign(order=df.set_index("platform")["fraction"])
    .reset_index(drop=False)
)
c2 = (
    alt.Chart(df)
    .mark_bar()
    .encode(
        y=alt.Y(
            "cell_type_coarse",
            title=None,
            sort=df.loc[lambda x: x["platform"] == "10x"]
            .sort_values("fraction", ascending=False)["cell_type_coarse"]
            .tolist(),
        ),
        x=alt.X(
            "fraction",
            sort=df.sort_values("fraction")["platform"].tolist(),
            scale=alt.Scale(domain=[0, 1]),
        ),
        color=alt.Color(
            "platform",
            scale=sh.colors.altair_scale("platform"),
            sort=alt.EncodingSortField("order"),
        ),
        order=alt.Order("order", sort="descending"),
    )
) + c_vline

c1 & c2

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
def fmt_label(total_series, min_value=0):
    def fmt(pct):
        n = int(pct / 100 * sum(total_series))
        if n > min_value:
            return n
        else:
            return ""

    return fmt


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
    autopct=fmt_label(df["n_samples"], 7),
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
    ).assign(
        condition = lambda x: [{"NSCLC": "NOS", "healthy_control": "non-cancer"}.get(_, _) for _ in x["condition"]]
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
    autopct=fmt_label(df["n_patients"], 7),
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


fig, ax = plt.subplots()
ax.pie(
    df["n_datasets"],
    labels=df["platform"],
    autopct=fmt_label(df["n_datasets"], 0),
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

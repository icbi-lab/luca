# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

# %%
# %load_ext autoreload
# %autoreload 2
import scanpy as sc
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import scvi
import pandas as pd
import scanpy_helpers.de as de
import altair as alt
from IPython.core.display import display, HTML

sc.set_figure_params(figsize=(5, 5))

# %%
gene_sets = {
    "C-type Lectins": (
        """
        CLEC7A
        CLEC1A
        CLEC14A
        CLEC4E
        CLEC6A
        CLEC4C
        CLEC5A
        CD209
        CLEC4A
        CLEC12A
        CLEC9A
        CLEC12B
        OPRL1
        CLEC3B
        CLEC4A
        CLEC2D
        KLRB1
        MRC1
        KLRD1
        WDFY4
        """
    ).split(),
    "TLRs": (
        """
        TLR1
        TLR2
        TLR3
        TLR4
        TLR5
        TLR6
        TLR7
        TLR8
        TLR10
        """
    ).split(),
    "Intracellulaar PRR": (
        """
        STING1
        CGAS
        ADA
        DDX58
        NOD1
        NOD2
        NLRP3
        NLRC4
        NLRP1
        """
    ).split(),
    "Downstream": (
        """
        TICAM1
        TICAM2
        AGER
        NFATC1
        NFATC2
        NFATC3
        NFATC4
        NFAT5
        TIRAP
        MYD88
        IRAK4
        IRAK1
        TRAF3
        TRAF6
        IKBKE
        IKBKB
        MAP3K7
        IRF7
        NFKB1
        NFKB2
        RIPK2
        ICE1
        IFIH1
        IRF3
        DHX58
        NWD1
        """
    ).split(),
    "Downstream Cytokines/ILs": (
        """
        CCL2
        CCL5
        CCL7
        CXCL8
        CXCL10
        CXCR6
        IL6
        IL12B
        IL23R
        IL17RA
        IL1B
        IL18
        IFNG
        IL2
        PLD3
        PLD4
        IL5
        IL10
        """
    ).split(),
    "STING-related genes": (
        """
        STING1
        CGAS
        IKBKB
        IRF3
        STAT6
        TBK1
        TRAF6
        TRAF3
        NFKB1
        NFKB2
        IFNB1
        CCL7
        CCL2
        CCL20
        DDX41
        IFI16
        ISG15
        ISG20
        ISG20L2
        IL6
        IL7
        IL1A
        IL1B
        CXCL1
        CXCL10
        CXCL8
        IFIT1
        PLD3
        """
    ).split(),
}

# %%
res_dir = "../../data/70_plots/71_de_analysis"

# %%
adata = sc.read_h5ad("../../data/60_infercnv/all-annotated-integrated.h5ad")
adata_scvi = sc.read_h5ad(
    "../../data/50_integrate_scrnaseq_data/52_run_scanvi/integrated_merged_all.h5ad"
)

# %%
# %%capture
scvi_model = scvi.model.SCANVI.load(
    "../../data/50_integrate_scrnaseq_data/52_run_scanvi/scvi_model_merged_all/",
    adata=adata_scvi,
)

# %%
adata.obs["condition"].unique()

# %%
adata.obs["disease_state"] = "other"
adata.obs.loc[
    adata.obs["condition"] == "healthy_control", "disease_state"
] = "healthy_control"
adata.obs.loc[
    adata.obs["origin"].isin(["tumor_primary", "tumor_metastasis"]), "disease_state"
] = "tumor"
adata.obs.loc[adata.obs["condition"] == "COPD", "disease_state"] = "COPD"

# %%
adata.obs.groupby(["dataset", "disease_state"]).size().reset_index(name="count")

# %%
adata.obs["disease_state"].value_counts()

# %%
adata = adata[(adata.obs["disease_state"] != "other"), :]

# %%
# Cell-types of interest
cell_types_coarse = {
    "B cell": True,
    "cDC1": True,
    "cDC2": True,
    "Ciliated": True,
    "mDC mature": True,
    "Endothelial cell": True,
    "Endothelial cell lymphatic": True,
    "Epithelial cell": [
        "Epithelial cell other (benign)",
        "Club/Goblet",
        "Alveolar cell type 1",
        "Alveolar cell type 2",
    ],
    "Fibroblast": ["Fibroblast adventitial", "Fibroblast", "Fibroblast alveolar"],
    "stromal (non-Fibroblast)": [
        "Smooth muscle cell",
        "Pericyte",
        "Mesothelial",
        "stromal other",
    ],
    "Monocyte": ["Monocyte non-conventional", "Monocyte conventional"],
    "Macrophage": ["Macrophage", "Macrophage FABP4+"],
    "Plasma cell": True,
    "Mast cell": ["other"],
    "myeloid (other)": ["myeloid other"],
    "T cell": ["T cell CD8", "T cell CD4", "T cell dividing", "T reg"],
    "NK cell": True,
    "Epithelial cell (malignant)": True,
}

# %%
adata.obs["cell_type_coarse"] = "other"
for new_label, cell_types in cell_types_coarse.items():
    if cell_types is True:
        cell_types = [new_label]
    for tmp_cell_type in cell_types:
        mask = adata.obs["cell_type"] == tmp_cell_type
        assert np.sum(mask), f"{tmp_cell_type} not found!"
        adata.obs.loc[mask, "cell_type_coarse"] = new_label

# %%
sc.pl.umap(adata, color="cell_type_coarse")

# %%
adata.obs["cell_type_disease_state"] = [
    f"{ct}_{ds}"
    for ct, ds in zip(adata.obs["cell_type_coarse"], adata.obs["disease_state"])
]

# %%
adata_scvi.obs["disease_state"] = adata.obs["disease_state"]
adata_scvi.obs["cell_type_coarse"] = adata.obs["cell_type_coarse"]
adata_scvi.obs["cell_type_disease_state"] = adata.obs["cell_type_disease_state"]
adata_scvi.obs["cell_type_dataset_disease_state"] = [
    f"{cell_type}_{dataset}_{disease_state}"
    for cell_type, dataset, disease_state in zip(
        adata_scvi.obs["cell_type_coarse"],
        adata_scvi.obs["dataset"],
        adata_scvi.obs["disease_state"],
    )
]

# %%
results = []
for cell_type in adata.obs["cell_type_coarse"].unique():
    res = scvi_model.differential_expression(
        adata_scvi[adata_scvi.obs["cell_type_coarse"] == cell_type, :],
        groupby="disease_state",
        group1=["tumor", "COPD"],
        group2="healthy_control",
        batch_correction=True,
        mode="change",
    ).assign(cell_type=cell_type)
    results.append(res)

# %%
res_all = pd.concat(results)

# %%
res_all.to_csv(f"{res_dir}/de_scvi.csv")

# %%
res_all.loc[
    (res_all["raw_normalized_mean1"] > 0.5) | (res_all["raw_normalized_mean2"] > 0.5), :
].to_csv(f"{res_dir}/de_scvi_filtered.csv")

# %% [markdown]
# ## pseudobulk-Analyses

# %%
ad = de.make_pseudobulk(
    adata_scvi,
    groupby="cell_type_dataset_disease_state",
    include_variables=["cell_type_coarse", "dataset", "disease_state"],
)

# %%
ad.obs["cell_threshold"] = ad.obs["n_cells"] >= 100

# %%
sc.pp.normalize_total(ad, target_sum=1e6)


# %%
def make_chart(ad, cell_type, gene_set):
    tmp_adata = ad[ad.obs["cell_type_coarse"] == cell_type, gene_sets[gene_set]]
    tmp_df = (
        pd.DataFrame(
            tmp_adata.X, index=tmp_adata.obs_names, columns=tmp_adata.var_names
        )
        .join(tmp_adata.obs)
        .melt(id_vars=tmp_adata.obs.columns, var_name="gene_symbol", value_name="cpm")
    )
    return (
        (
            alt.Chart(tmp_df)
            .mark_boxplot()
            .encode(x="disease_state", y="cpm")
            .properties(height=150)
            + (
                alt.Chart(tmp_df)
                .mark_circle(size=60, fillOpacity=1)
                .encode(
                    alt.X("disease_state", title=None),
                    y="cpm",
                    color=alt.Color("dataset", scale=alt.Scale(scheme="paired")),
                    stroke=alt.Stroke(
                        "cell_threshold",
                        scale=alt.Scale(
                            domain=["true", "false"], range=["black", "lightgrey"]
                        ),
                    ),
                )
            )
        )
        .facet("gene_symbol")
        .resolve_scale(y="independent")
        .properties(title=f"{cell_type}: {gene_set}")
    )


# %%
for gene_set in gene_sets.keys():
    display(HTML(f"<h3>{gene_set}</h3>"))
    display(
        alt.vconcat(
            *[
                make_chart(ad, ct, gene_set)
                for ct in ad.obs["cell_type_coarse"].unique()
                if not pd.isnull(ct)
            ]
        ).configure_mark(opacity=0.9)
    )

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
import pandas as pd
import scanpy_helpers as sh
from scanpy_helpers.compare_groups.pl import plot_lm_result_altair
import scanpy as sc
import numpy as np
from nxfvars import nxfvars

# %%
path_prefix = nxfvars.get(
    "path_prefix",
    "../../data/30_downstream_analyses/plots_and_comparisons/91_compare_groups/{comparison}/artifacts",
)
deseq2_path_prefix = nxfvars.get(
    "deseq2_path_prefix",
    "../../data/30_downstream_analyses/de_analysis/{comparison}/de_deseq2",
)

# %%
tools = ["dorothea", "progeny", "cytosig"]
comparisons = [
    "patient_infiltration_status",
    "patient_immune_infiltration",
    "patient_immune_infiltration_treatment_coding",
    "patient_immune_infiltration_treatment_coding_random",
    "luad_lusc",
    "early_advanced",
]

# %%
results = {}
for comparison in comparisons:
    results[comparison] = {}
    for tool in tools:
        try:
            results[comparison][tool] = pd.read_csv(
                (path_prefix + "/{comparison}_{tool}.tsv").format(
                    comparison=comparison, tool=tool
                ),
                sep="\t",
            )
        except IOError:
            raise

# %% [markdown]
# # LUAD vs LUSC

# %% [markdown]
# ## Dorothea
# ### Tumor cells

# %%
results["luad_lusc"]["dorothea"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="TFs (tumor cells LUAD/LUSC)")

# %% [markdown]
# ## Progeny
# ### Tumor cells

# %%
results["luad_lusc"]["progeny"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1)

# %% [markdown]
# ## Cytosig
# ### Tumor cells

# %%
results["luad_lusc"]["cytosig"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Cytosig (tumor cells)")

# %% [markdown]
# ### Stromal cells

# %%
results["luad_lusc"]["cytosig"].loc[lambda x: x["cell_type"] == "Stromal", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Cytosig (stromal cells)")

# %% [markdown]
# # Early advanced

# %% [markdown]
# ## Progeny
# ### Tumor cells

# %%
results["early_advanced"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1
)

# %% [markdown]
# ## Cytosig
# (Tumor and stromal cells)

# %%
results["early_advanced"]["cytosig"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Cytosig (Neutrophils cells)"
)

# %%
results["early_advanced"]["cytosig"].loc[lambda x: x["cell_type"] == "Stromal", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Cytosig (Neutrophils cells)")

# %% [markdown] tags=[]
# # Infiltration groups

# %%
de_res_tumor_cells = (
    pd.concat(
        [
            pd.read_csv(
                (
                    deseq2_path_prefix
                    + "/{comparison}_primary_tumor_adata_primary_tumor_tumor_cells_DESeq2_result.tsv"
                ).format(comparison=k),
                sep="\t",
            ).assign(group=k.upper())
            for k in ["t_desert", "m_desert", "b_desert"]
        ]
    )
    .fillna(1)
    .pipe(sh.util.fdr_correction)
    .drop("gene_id", axis="columns")
    .rename(columns={"gene_id.1": "gene_id"})
)

# %%
de_res_tumor_cells

# %%
selected_genes = (
    de_res_tumor_cells.groupby("group").apply(lambda x: x.head(15))["gene_id"].values
)

# %%
de_res_tumor_cells.loc[lambda x: x["gene_id"].isin(selected_genes)].assign(
    log2FoldChange=lambda x: np.clip(x["log2FoldChange"], -5, 5)
).pipe(
    plot_lm_result_altair,
    title="Differential genes (tumor cells)",
    color="log2FoldChange",
    x="gene_id",
)

# %% [markdown]
# ## Dorothea

# %%
results["patient_immune_infiltration_treatment_coding"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair,
    title="Differential TFs (tumor cells)",
)

# %% [markdown]
# ## Progeny

# %%
results["patient_immune_infiltration_treatment_coding"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1
)

# %% [markdown]
# ## Cytosig
# ### Infiltration subtypes

# %%
tmp_cytosig = (
    results["patient_immune_infiltration_treatment_coding"]["cytosig"]
    .loc[lambda x: x["cell_type"].isin(["Tumor cells"]), :]
    .pipe(sh.util.fdr_correction)
)

# %%
for ct in tmp_cytosig["cell_type"].unique():
    try:
        plot_lm_result_altair(
            tmp_cytosig.loc[lambda x: x["cell_type"] == ct],
            title=f"Cytosig for {ct}",
            p_cutoff=0.1,
        ).display()
    except AttributeError:
        pass

# %% [markdown] tags=[]
# # Infiltration groups (sum2zero coding)

# %% [markdown]
# ## Dorothea

# %%
results["patient_immune_infiltration"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair,
    title="Differential TFs (tumor cells)",
)

# %% [markdown]
# ## Progeny

# %%
results["patient_immune_infiltration"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1
)

# %% [markdown]
# ## Cytosig
# ### Infiltration subtypes

# %%
tmp_cytosig = (
    results["patient_immune_infiltration"]["cytosig"]
    .loc[lambda x: x["cell_type"].isin(["Tumor cells"]), :]
    .pipe(sh.util.fdr_correction)
)

# %%
for ct in tmp_cytosig["cell_type"].unique():
    try:
        plot_lm_result_altair(
            tmp_cytosig.loc[lambda x: x["cell_type"] == ct],
            title=f"Cytosig for {ct}",
            p_cutoff=0.1,
        ).display()
    except AttributeError:
        pass

# %% [markdown]
# ---

# %% [markdown]
# # Infiltration groups (random control)

# %% [markdown]
# ## Dorothea

# %%
results["patient_immune_infiltration_treatment_coding_random"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair,
    title="Differential TFs (tumor cells)",
)

# %% [markdown]
# ## Progeny

# %%
results["patient_immune_infiltration_treatment_coding_random"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1
)

# %% [markdown]
# ## Cytosig
# ### Infiltration subtypes

# %%
tmp_cytosig = (
    results["patient_immune_infiltration_treatment_coding_random"]["cytosig"]
    .loc[lambda x: x["cell_type"].isin(["Tumor cells"]), :]
    .pipe(sh.util.fdr_correction)
)

# %%
for ct in tmp_cytosig["cell_type"].unique():
    try:
        plot_lm_result_altair(
            tmp_cytosig.loc[lambda x: x["cell_type"] == ct],
            title=f"Cytosig for {ct}",
            p_cutoff=0.1,
        ).display()
    except AttributeError:
        pass

# %%

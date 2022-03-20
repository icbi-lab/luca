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

# %%
tools = ["dorothea", "progeny", "cytosig", "cpdb"]
comparisons = [
    "tumor_normal",
    "infiltration_status",
    "infiltration_type",
    "patient_immune_infiltration",
    "patient_immune_infiltration_condition",
    "patient_immune_infiltration_treatment_coding",
    "patient_immune_infiltration_treatment_coding_condition",
    "patient_immune_infiltration_treatment_coding_condition_random",
    "luad_lusc",
    "early_advanced",
    "early_advanced_condition",
]

# %%
results = {}
for comparison in comparisons:
    results[comparison] = {}
    for tool in tools:
        try:
            results[comparison][tool] = pd.read_csv(
                f"../../data/30_downstream_analyses/plots_and_comparisons/91_compare_groups/artifacts/{comparison}_{tool}.tsv",
                sep="\t",
            )
        except IOError:
            pass

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

# %%
results["early_advanced_condition"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1
)

# %% [markdown]
# ## Cytosig
# (Tumor and stromal cells)

# %%
results["early_advanced_condition"]["cytosig"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Cytosig (Neutrophils cells)"
)

# %%
results["early_advanced_condition"]["cytosig"].loc[
    lambda x: x["cell_type"] == "Stromal", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Cytosig (Neutrophils cells)"
)

# %% [markdown]
# # Infiltration groups

# %% [markdown]
# ## Dorothea

# %%
results["patient_immune_infiltration_treatment_coding_condition"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair,
    title="Differential TFs (tumor cells)",
)

# %% [markdown]
# ## Progeny

# %%
results["patient_immune_infiltration_treatment_coding_condition"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1
)

# %% [markdown]
# ## Cytosig
# ### Infiltration subtypes

# %%
tmp_cytosig = (
    results["patient_immune_infiltration_treatment_coding_condition"]["cytosig"]
    .loc[
        lambda x: x["cell_type"].isin(["Tumor cells"]), :
    ]
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
results["patient_immune_infiltration_treatment_coding_condition_random"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair,
    title="Differential TFs (tumor cells)",
)

# %% [markdown]
# ## Progeny

# %%
results["patient_immune_infiltration_treatment_coding_condition_random"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1
)

# %% [markdown]
# ## Cytosig
# ### Infiltration subtypes

# %%
tmp_cytosig = (
    results["patient_immune_infiltration_treatment_coding_condition_random"]["cytosig"]
    .loc[
        lambda x: x["cell_type"].isin(["Tumor cells"]), :
    ]
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

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
    "patient_immune_infiltration_treatment_coding_condition2",
    "luad_lscc",
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
# ## Dorothea

# %% [markdown]
# ### Tumor cells

# %%
results["luad_lscc"]["dorothea"].loc[lambda x: x["cell_type"] == "Neutrophils", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="TFs (tumor cells LUAD/LSCC)")

# %%
results["luad_lscc"]["dorothea"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="TFs (tumor cells LUAD/LSCC)")

# %% [markdown]
# ### Neutrophils

# %%
results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].to_csv("/home/sturm/Downloads/neutrophils_tfs.tsv", sep="\t")

# %%
results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="TFs (Neutrophils tumor/normal)"
)

# %%
results["early_advanced"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="TFs (Neutrophils tumor/normal)"
)

# %% [markdown]
# ## Progeny

# %%
results["patient_immune_infiltration_treatment_coding"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)"
)

# %%
results["patient_immune_infiltration_treatment_coding_condition"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)"
)

# %%
results["patient_immune_infiltration_treatment_coding_condition2"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Differential pathways (tumor cells)"
)

# %%
results["luad_lscc"]["progeny"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Differential pathways (tumor cells)", p_cutoff=1)

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
# ### Infiltration subtypes

# %%
tmp_cytosig = (
    results["patient_immune_infiltration_treatment_coding"]["cytosig"]
    .loc[
        lambda x: x["cell_type"].isin(["Tumor cells", "Stromal", "Endothelial cell"]), :
    ]
    .pipe(sh.util.fdr_correction)
)

# %%
for ct in tmp_cytosig["cell_type"].unique():
    try:
        plot_lm_result_altair(
            tmp_cytosig.loc[lambda x: x["cell_type"] == ct], title=f"Cytosig for {ct}"
        ).display()
    except AttributeError:
        pass

# %%
tmp_cytosig = (
    results["patient_immune_infiltration_treatment_coding_condition"]["cytosig"]
    .loc[
        lambda x: x["cell_type"].isin(["Tumor cells", "Stromal", "Endothelial cell"]), :
    ]
    .pipe(sh.util.fdr_correction)
)

# %%
for ct in tmp_cytosig["cell_type"].unique():
    try:
        plot_lm_result_altair(
            tmp_cytosig.loc[lambda x: x["cell_type"] == ct], title=f"Cytosig for {ct}"
        ).display()
    except AttributeError:
        pass

# %%
tmp_cytosig = (
    results["patient_immune_infiltration_treatment_coding_condition2"]["cytosig"]
    .loc[
        lambda x: x["cell_type"].isin(["Tumor cells", "Stromal", "Endothelial cell"]), :
    ]
    .pipe(sh.util.fdr_correction)
)

# %%
for ct in tmp_cytosig["cell_type"].unique():
    try:
        plot_lm_result_altair(
            tmp_cytosig.loc[lambda x: x["cell_type"] == ct], title=f"Cytosig for {ct}"
        ).display()
    except AttributeError:
        pass

# %% [markdown]
# ### LUAD / LUSC

# %%
ad_cyto_luad_lusc = sc.read_h5ad("../../data/30_downstream_analyses/plots_and_comparisons/91_compare_groups/artifacts/luad_lscc_cytosig/Tumor cells.h5ad")

# %%
pb_cyto_laud_lusc = sh.pseudobulk.pseudobulk(ad_cyto_luad_lusc, groupby=["patient", "condition", "dataset"])

# %%
sh.pairwise.plot_paired(pb_cyto_laud_lusc, "condition", hue="dataset", var_names=["MCSF", "BMP2", "BMP4"], size=2, panel_size=(4, 4))

# %%
results["luad_lscc"]["cytosig"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Cytosig (tumor cells)")

# %%
results["luad_lscc"]["cytosig"].loc[lambda x: x["cell_type"] == "Stromal", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Cytosig (stromal cells)")

# %%
results["luad_lscc"]["cytosig"].loc[lambda x: x["cell_type"] == "Neutrophils", :].pipe(
    sh.util.fdr_correction
).pipe(plot_lm_result_altair, title="Cytosig (Neutrophils cells)")

# %% [markdown]
# ### Early/Advanced

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

# %%
results["early_advanced_condition"]["cytosig"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].pipe(sh.util.fdr_correction).pipe(
    plot_lm_result_altair, title="Cytosig (Neutrophils cells)"
)

# %%

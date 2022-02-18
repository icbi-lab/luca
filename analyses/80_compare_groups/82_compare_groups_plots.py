# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Differentially expressed dorothea TFs
#
#  * nothing significant for Neutrophils

# %%
results["luad_lscc"]["dorothea"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(sh.lm.plot_lm_result_altair, title="TFs (tumor cells LUAD/LSCC)")

# %%
results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].to_csv("/home/sturm/Downloads/neutrophils_tfs.tsv", sep="\t")

# %%
results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
].pipe(sh.util.fdr_correction).pipe(
    sh.lm.plot_lm_result_altair, title="TFs (Neutrophils tumor/normal)"
)

# %%
pb_dorothea = sh.pseudobulk.pseudobulk(
    datasets["tumor_normal"]["dorothea"]["Neutrophils"],
    groupby=["dataset", "patient", "origin"],
    aggr_fun=np.mean,
)

# %%
tfoi = results["tumor_normal"]["dorothea"].loc[
    lambda x: x["cell_type"] == "Neutrophils", :
]["variable"][:30]

# %%
sh.pairwise.plot_paired_fc(
    pb_dorothea, groupby="origin", paired_by="patient", metric="diff", var_names=tfoi
).properties(height=150)

# %%
sh.pairwise.plot_paired(
    pb_dorothea, groupby="origin", paired_by="patient", var_names=tfoi
)

# %% [markdown]
# ## Differentially expressed progeny pathways in tumor cells

# %%
results["patient_immune_infiltration_treatment_coding"]["progeny"].loc[
    lambda x: x["cell_type"] == "Tumor cells", :
].pipe(sh.util.fdr_correction).pipe(
    sh.lm.plot_lm_result_altair, title="Differential pathways (tumor cells)"
)

# %%
results["luad_lscc"]["progeny"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(sh.lm.plot_lm_result_altair, title="Differential pathways (tumor cells)")

# %% [markdown]
# ## Differential cytokine signalling in selected cell-types

# %%
results["patient_immune_infiltration_treatment_coding"]["cytosig"]["group"].unique()

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
        sh.lm.plot_lm_result_altair(
            tmp_cytosig.loc[lambda x: x["cell_type"] == ct], title=f"Cytosig for {ct}"
        ).display()
    except AttributeError:
        pass

# %%
results["luad_lscc"]["cytosig"].loc[lambda x: x["cell_type"] == "Tumor cells", :].pipe(
    sh.util.fdr_correction
).pipe(sh.lm.plot_lm_result_altair, title="Cytosig (tumor cells)")

# %% [markdown]
# ---

# %%
marker_genes = {}
for tf in regulons.columns:
    marker_genes[tf] = regulons[tf][lambda x: x != 0]

# %%
pd.concat(marker_genes).reset_index(name="direction").rename(
    columns={"level_0": "TF"}
).to_csv("/home/sturm/Downloads/dorothea_signature.csv")

# %%
marker_genes["SOX2"].index.tolist()

# %%
regulons["SOX2"][lambda x: x != 0]

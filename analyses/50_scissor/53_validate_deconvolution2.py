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
import scanpy_helpers as sh
from nxfvars import nxfvars
import pandas as pd
import scanpy as sc
from scanpy_helpers.compare_groups.pl import plot_lm_result_altair
import numpy as np
from lifelines import CoxPHFitter

# %%
immunedeconv_results_path = nxfvars.get(
    "immunedeconv_results",
    "../../data/13_tcga/immunedeconv/nsclc_deconvolution_remapped.tsv",
)
clinical_data_path = nxfvars.get(
    "clincal_data", "../../tables/tcga/clinical_data_for_scissor.tsv"
)

# %% [markdown]
# # Validate scissor results with deconvolution
#
# Pre-computed immunedeconv results are loaded from a TSV file. 

# %%
immunedeconv_results = pd.read_csv(immunedeconv_results_path, sep="\t").loc[
    lambda x: x["method"] != "cibersort", :
]
clinical_data = pd.read_csv(
    clinical_data_path, sep="\t", index_col="TCGA_patient_barcode"
)

# %%
clinical_data

# %%
immunedeconv_results

# %%
# Normalize each cell-type score between 0 and 1 such they are comparable between methods (and within the same method)
# We are looking for relative changes.
immunedeconv_results_norm = (
    (immunedeconv_results.dropna().set_index(["cell_type", "method"]))
    .apply(lambda row: row / np.max(row), axis=1)
    .reset_index()
)

# %%
# Convert to AnnData such that we can use scanpy helper functions for plotting
# and fittin the LM
adatas = {}
for method in immunedeconv_results["method"].unique():
    tmp_adata = sc.AnnData(
        immunedeconv_results_norm.loc[lambda x: x["method"] == method, :]
        .drop("method", axis="columns")
        .set_index("cell_type")
        .T
    )
    tmp_adata.obs = clinical_data.loc[tmp_adata.obs_names, :]
    adatas[method] = tmp_adata

# %% [markdown]
# ## LUAD vs LUSC

# %%
tmp_results = []
for method, tmp_adata in adatas.items():
    tmp_results.append(
        sh.compare_groups.lm.test_lm(
            tmp_adata,
            groupby="type",
            formula="~ C(type, Treatment('LUAD')) + tumor_stage_ajcc + gender + age",
            contrasts="Treatment('LUAD')",
            progress=False,
        ).assign(method=method)
    )
res_luad_lusc = pd.concat(tmp_results).pipe(sh.util.fdr_correction)

# %%
plot_lm_result_altair(
    res_luad_lusc, y="method", color="coef", title="LUAD vs LUSC (red = up in LUSC)"
)

# %% [markdown]
# ## early vs advanced

# %%
tmp_results = []
for method, tmp_adata in adatas.items():
    tmp_results.append(
        sh.compare_groups.lm.test_lm(
            tmp_adata,
            groupby="tumor_stage",
            formula="~ C(tumor_stage, Treatment(0.0)) + type + gender + age",
            contrasts="Treatment(0.0)",
            progress=False,
        ).assign(method=method)
    )
res_early_advanced = pd.concat(tmp_results).pipe(sh.util.fdr_correction)

# %%
plot_lm_result_altair(
    res_early_advanced,
    y="method",
    color="coef",
    p_cutoff=1,
    title="early vs advanced (red = up in advanced)",
)

# %% [markdown]
# ## Mutations

# %%
columns = ["kras_mutation", "egfr_mutation", "stk11_mutation", "tp53_mutation"]


# %%
def plot_mutation_luad_lusc(col):
    """Plot mutation association for LUAD and LUSC separately"""
    for tumor_type in ["LUAD", "LUSC"]:
        tmp_results = []
        for method, tmp_adata in adatas.items():
            tmp_adata2 = tmp_adata[tmp_adata.obs["type"] == tumor_type, :].copy()
            tmp_results.append(
                sh.compare_groups.lm.test_lm(
                    tmp_adata2,
                    groupby=col,
                    formula=f"~ C({col}, Treatment(1)) + tumor_stage_ajcc + gender + age",
                    contrasts="Treatment(1)",
                    progress=False,
                ).assign(method=method)
            )
        res = (
            pd.concat(tmp_results)
            .assign(pvalue=lambda x: x["pvalue"].fillna(1))
            .pipe(sh.util.fdr_correction)
        )
        plot_lm_result_altair(
            res,
            y="method",
            color="coef",
            p_cutoff=1,
            title=f"{col}: {tumor_type} (blue = up in mutated)",
        ).display()


# %%
for col in columns:
    plot_mutation_luad_lusc(col)

# %% [markdown]
# ## Survival

# %%
res_luad_lusc.head()


# %%
def test_coxph(
    adata,
    *,
    formula="",
    duration_col="time",
    event_col="status",
    strata=None,
    max_days=365 * 10,
):
    """
    Fit a coxph model

    Parameters
    ----------
    formula
        covariates to add to the formula. Can be an empty string, otherwise needs to start with "+"
    """
    df = adata.obs.join(
        pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    )
    res_builder = []
    for var in adata.var_names:
        tmp_df = df.copy()

        # subset to used columns and drop NAs
        use_cols = [x.strip() for x in formula.split("+")]
        use_cols = (
            set(use_cols)
            | (set(strata) if strata is not None else set())
            | {duration_col, event_col, var}
        ) - {""}
        tmp_df = tmp_df.loc[:, use_cols].dropna()
        tmp_df = tmp_df.loc[lambda x: x[duration_col] <= max_days, :]

        # sanitize variable name because the formula doesn't support quoting (or I don't know how)
        var_san = var.replace(" ", "_").replace("+", "")
        tmp_df = tmp_df.rename(columns={var: var_san})

        # fit the coxph model
        cph = CoxPHFitter()
        cph.fit(
            tmp_df,
            duration_col=duration_col,
            event_col=event_col,
            formula=f"{var_san} {formula}",
            strata=strata,
        )

        # extract relevant results
        tmp_res = cph.summary.loc[var_san, ["coef", "p"]].to_dict()
        tmp_res["variable"] = var
        res_builder.append(tmp_res)

    return pd.DataFrame(res_builder).rename(columns={"p": "pvalue"})


# %% [markdown]
# ### with covariates

# %%
for tumor_type in ["all", "LUAD", "LUSC"]:
    tmp_results = []
    for method, tmp_adata in adatas.items():
        if tumor_type != "all":
            tmp_adata = tmp_adata[tmp_adata.obs["type"] == tumor_type, :].copy()
        tmp_results.append(
            test_coxph(
                tmp_adata,
                formula="+ age",
                strata=["tumor_stage_ajcc", "gender"],
            ).assign(method=method)
        )
    res_survival = pd.concat(tmp_results).pipe(sh.util.fdr_correction)
    plot_lm_result_altair(
        res_survival,
        y="method",
        color="coef",
        p_cutoff=1,
        title=f"survival: {tumor_type} (red = worse survival)",
        value_max=3,
    ).display()

# %% [markdown]
# ### without covariates

# %%
for tumor_type in ["all", "LUAD", "LUSC"]:
    tmp_results = []
    for method, tmp_adata in adatas.items():
        if tumor_type != "all":
            tmp_adata = tmp_adata[tmp_adata.obs["type"] == tumor_type, :].copy()
        tmp_results.append(
            test_coxph(
                tmp_adata,
                # formula="+ age",
                # strata=["tumor_stage_ajcc", "gender"],
            ).assign(method=method)
        )
    res_survival = pd.concat(tmp_results).pipe(sh.util.fdr_correction)
    plot_lm_result_altair(
        res_survival,
        y="method",
        color="coef",
        p_cutoff=1,
        title=f"survival: {tumor_type} (red = worse survival)",
        value_max=3,
    ).display()

# %%

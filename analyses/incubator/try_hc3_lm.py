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
import pandas as pd
from tqdm import tqdm
import numpy as np
import statsmodels.formula.api as smf

# %%
adata_cpdb = sc.read_h5ad(
    "../../data/30_downstream_analyses/cell2cell_major/cpdb_h5ad/artifacts/adata_cpdb.h5ad"
)

# %%
adata_cpdb.obs["immune_infiltration"] = (
    patient_stratification.loc[:, ["patient", "immune_infiltration"]]
    .merge(adata_cpdb.obs.reset_index())
    .set_index("index")
    .loc[:, "immune_infiltration"]
)

# %%
adata = sc.read_h5ad(
    "../../data/30_downstream_analyses/03_update_annotation/artifacts/full_atlas_merged.h5ad"
)

# %%
adata_cpdb

# %%
adata_cpdb.var

# %%
plasma = adata_cpdb[
    (adata_cpdb.obs["origin"] == "tumor_primary")
    & (~adata_cpdb.obs["immune_infiltration"].isnull()),
    (adata_cpdb.var["target"] == "KIR2DL3")
    & (adata_cpdb.var["cluster_2"] == "Plasma cell"),
]

# %%
plasma.obs

# %%
sc.pl.dotplot(
    adata[adata.obs["origin"] == "tumor_primary", :],
    var_names=["KIR2DL3"],
    groupby=["cell_type_coarse"],
)

# %%
sc.pl.dotplot(
    adata[
        (adata.obs["origin"] == "tumor_primary")
        & (adata.obs["cell_type_coarse"] == "Plasma cell"),
        :,
    ],
    var_names=["KIR2DL3"],
    groupby=["patient"],
)

# %%
sc.pl.dotplot(plasma, var_names=plasma.var_names, groupby="patient")

# %%
sc.pl.dotplot(plasma, var_names=plasma.var_names, groupby="immune_infiltration")

# %%
sc.pl.heatmap(plasma, var_names=plasma.var_names, groupby="immune_infiltration")


# %%
def get_stat(ad, b_mask, desert_mask):
    return np.mean(ad[b_mask].X, axis=0) - np.mean(ad[desert_mask].X, axis=0)


# %%
tmp_plasma = plasma.copy() # plasma[plasma.obs["immune_infiltration"].isin(["B", "desert"]), :].copy()
stats = []
group_vec = tmp_plasma.obs["immune_infiltration"].values.copy()
for i in tqdm(range(1000)):
    np.random.shuffle(group_vec)
    stats.append(get_stat(tmp_plasma, group_vec == "B", group_vec == "desert"))

# %%
stat = get_stat(tmp_plasma, tmp_plasma.obs["immune_infiltration"] == "B", tmp_plasma.obs["immune_infiltration"] == "desert")

# %%
stat

# %%
(np.sum((np.vstack(stats) > stat), axis=0) + 1) / (len(stats) + 1)

# %%
df = plasma.obs.join(pd.DataFrame(plasma.X, index=plasma.obs_names, columns=plasma.var_names))

# %%
df

# %%
var = "HLA-C_KIR2DL3_Alveolar cell type 2_Plasma cell"
mod = smf.ols(f"Q('{var}') ~ C(immune_infiltration, Sum) + dataset + tumor_stage", df)

# %%
res = mod.fit()

# %%
sm = res.summary()

# %%
res_hc3 = res.get_robustcov_results(cov_type="HC3")

# %%
res.pvalues

# %%
type(res)

# %%
type(res_hc3)

# %%
import statsmodels

# %%
statsmodels.regression.linear_model.RegressionResultsWrapper(res_hc3).params

# %%
res_hc3.pvalues

# %%
res_hc3.params

# %%
res_hc3.params

# %%
res_hc3.summary()

# %%
res_hc3.params

# %%
res.params

# %%
res_rlm = smf.rlm(f"Q('{var}') ~ C(immune_infiltration, Sum) + dataset + tumor_stage", df).fit()

# %%
res.summary()

# %%
res_hc3 = res.get_robustcov_results(cov_type="HC3")

# %%
res_hc3.summary()

# %%
res_hc3.params

# %% [markdown]
# ----

# %%
nk = adata_cpdb[
    (adata_cpdb.obs["origin"] == "tumor_primary")
    & (~adata_cpdb.obs["immune_infiltration"].isnull()),
    (adata_cpdb.var["target"] == "KIR2DL3")
    & (adata_cpdb.var["cluster_2"] == "NK cell"),
]

# %%

# %%
sc.pl.dotplot(nk, var_names=nk.var_names, groupby="patient")

# %%
sc.pl.dotplot(nk, var_names=nk.var_names, groupby="immune_infiltration")

# %%
sc.pl.heatmap(nk, var_names=nk.var_names, groupby="immune_infiltration")

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
adata = sc.read_h5ad(
    "../../data/30_downstream_analyses/de_analysis/b_desert/adata_by_cell_type/adata_primary_tumor_plasma_cell.h5ad"
)

# %%
patient_stratification = pd.read_csv(
    "../../data/30_downstream_analyses/stratify_patients/artifacts/patient_stratification.csv"
)

# %%
adata.obs

# %%
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# %%
sc.tl.pca(adata)

# %%
sc.pp.neighbors(adata)

# %%
sc.tl.umap(adata)

# %%
sc.pl.umap(adata, color=["TM4SF1", "S100A10", "MZB1"])

# %%
adata.obs

# %%
sc.tl.rank_genes_groups(
    adata, groupby="immune_infiltration", groups=["case"], reference="control"
)

# %%
sc.pl.rank_genes_groups(
    adata,
)

# %%
sc.pl.dotplot(adata, groupby=["immune infiltration

# %%
b2 = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/31_cell_types_coarse/by_cell_type/adata_cell_type_coarse_t_cell.umap_leiden.h5ad"
)

# %%
sc.pl.umap(
    b2,
    color=["CD3E", "NCAM1", "IGKC", "MZB1", "IGLC3", "IGHG4", "IGLC2"],
    cmap="inferno",
    size=20,
)

# %%

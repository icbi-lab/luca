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
import threadpoolctl
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st
import matplotlib
import scanpy_helpers as sh
import warnings
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map
import numpy as np
import scipy.stats

warnings.filterwarnings("ignore")

# %%
artifact_dir = "/home/sturm/Downloads"

# %%
adata = sc.read_h5ad(
    "../../data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad"
)

# %%
adata.obs["origin2"] = ["normal" if "normal" in x else x for x in adata.obs["origin"]]

# %%
goi = ["CBLB", "NR2F6"]


# %%
def prep_pseudobulk(ct):
    adata_tumor_normal = sh.pseudobulk.pseudobulk(
        adata[
            (adata.obs["cell_type_major"] == ct)
            & adata.obs["origin2"].isin(["tumor_primary", "normal"]),
            :,
        ].copy(),
        groupby=["dataset", "patient", "origin2"],
        min_obs=10,
    )
    sc.pp.normalize_total(adata_tumor_normal, target_sum=1e6)
    sc.pp.log1p(adata_tumor_normal)
    return adata_tumor_normal


# %%
pseudobulks = {
    ct: prep_pseudobulk(ct)
    for ct in tqdm(sorted(adata.obs["cell_type_major"].unique()))
}
del pseudobulks["other"]


# %%
def co_expression(ct, pb):
    corr = {g: {} for g in goi}
    pvals = {g: {} for g in goi}
    mean_expr = {g: {} for g in goi}
    mean_expr_target = {g: {} for g in goi}

    for query in pb.var_names:
        for ref in goi:
            corr[ref][query], pvals[ref][query] = scipy.stats.pearsonr(
                pb[:, ref].X[:, 0], pb[:, query].X[:, 0]
            )
            mean_expr[ref][query] = np.mean(pb[:, ref].X[:, 0])
            mean_expr_target[ref][query] = np.mean(pb[:, query].X[:, 0])

    return (
        pd.DataFrame.from_dict(corr)
        .join(pd.DataFrame.from_dict(pvals), rsuffix="_p")
        .join(pd.DataFrame.from_dict(mean_expr), rsuffix="_mean_query")
        .join(pd.DataFrame.from_dict(mean_expr_target), rsuffix="_mean_target")
        .assign(cell_type=ct)
    )


# %%
co_expr = pd.concat(
    process_map(co_expression, pseudobulks.keys(), pseudobulks.values())
)

# %%
co_expr_cblb = (
    co_expr.loc[
        co_expr.index != "CBLB",
        lambda x: x.columns.str.contains("CBLB") | (x.columns == "cell_type"),
    ]
    .dropna(how="any")
    .sort_values("CBLB")
    .loc[lambda x: (x["CBLB_mean_query"] >= 0.5) & (x["CBLB_mean_target"] >= 0.5)]
    .pipe(sh.util.fdr_correction, pvalue_col="CBLB_p")
)

# %%
co_expr_cblb.to_csv(f"{artifact_dir}/co_expression_filtered_cblb.csv")

# %%
co_expr_cblb

# %%
co_expr_nr2f6 = (
    co_expr.loc[
        co_expr.index != "NR2F6",
        lambda x: x.columns.str.contains("NR2F6") | (x.columns == "cell_type"),
    ]
    .dropna(how="any")
    .sort_values("NR2F6")
    .loc[lambda x: (x["NR2F6_mean_query"] >= 0.5) & (x["NR2F6_mean_target"] >= 0.5)]
    .pipe(sh.util.fdr_correction, pvalue_col="NR2F6_p")
)

# %%
co_expr_nr2f6.to_csv(f"{artifact_dir}/co_expression_filtered_nr2f6.csv")

# %%
co_expr_nr2f6


# %%
def plot_for_cell_type(ct, pb):
    deseq2_res_tumor_normal = pd.read_csv(
        f"../../data/30_downstream_analyses/de_analysis/tumor_normal/de_deseq2/adata_tumor_normal_{ct.lower().replace(' ', '_').replace('+', '_')}_DESeq2_result.tsv",
        sep="\t",
        index_col=0,
    ).set_index("gene_id.1")
    pvalues = deseq2_res_tumor_normal.reindex(goi)["pvalue"].values
    sh.pairwise.plot_paired(
        pb,
        groupby="origin2",
        paired_by="patient",
        var_names=goi,
        pvalues=pvalues,
        pvalue_template="DESeq2 p={:.2f}",
        ylabel="log CPM",
        show_legend=False,
    )


# %% [markdown]
# ## Expression of CBLB and NF2

# %%
for ct, pb in pseudobulks.items():
    print(ct)
    plot_for_cell_type(ct, pb)

# %% [markdown]
# ### overall expression

# %%
bdata = sh.pseudobulk.pseudobulk(
    adata[
        adata.obs["origin2"].isin(["tumor_primary"]),
        :,
    ].copy(),
    groupby=["dataset", "patient", "cell_type_major"],
    min_obs=10,
)

# %%
sc.pl.matrixplot(bdata, var_names=goi, groupby="cell_type_major", swap_axes=True)

# %%

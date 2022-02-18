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
import scanpy as sc
import pandas as pd
from nxfvars import nxfvars
import pickle
from pathlib import Path
from tqdm import tqdm

# %%
adata_atlas = nxfvars.get(
    "adata_atlas",
    "../../data/30_downstream_analyses/02_integrate_into_atlas/artifacts/full_atlas_merged.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/")
squidpy_dir = nxfvars.get(
    "squidpy_dir", "../../data/30_downstream_analyses/cell2cell/squidpy/"
)

# %% [markdown]
# ## Build anndata object of cell2cell interactions
#
# Load the results obtained for each sample and 
# load them into an anndata file (each obs is a sample and each var a pair of receptor/ligand)

# %%
adata = sc.read_h5ad(adata_atlas)

# %%
cpdb_res = {}
for f in Path(squidpy_dir).glob("**/*.pkl"):
    sample = f.name.replace(Path(adata_atlas).stem + "_", "").replace(".pkl", "")
    with open(f, "rb") as fh:
        cpdb_res[sample] = pickle.load(fh)

# %%
dfs_melt = {}
for k in tqdm(cpdb_res):
    dfs_melt[k] = (
        cpdb_res[k]["means"]
        .reset_index()
        .melt(id_vars=["source", "target"], value_name=k)
    )

# %%
var = pd.concat(
    [
        df.loc[lambda x: x[k] != 0, ["source", "target", "cluster_1", "cluster_2"]]
        for k, df in tqdm(dfs_melt.items())
    ]
).drop_duplicates()

# %%
var = var.assign(idx=lambda x: ["_".join(t[1:]) for t in x.itertuples()]).set_index(
    "idx"
)

# %%
for k, df in tqdm(dfs_melt.items()):
    tmp_series = (
        df.loc[lambda x: x[k] != 0, :]
        .assign(idx=lambda x: ["_".join(t[1:-1]) for t in x.itertuples()])
        .set_index("idx")[k]
    )
    var[k] = tmp_series

# %%
ad_cpdb = sc.AnnData(var=var.iloc[:, :4], X=var.iloc[:, 4:].T.fillna(0))

# %%
sample_info = (
    adata.obs.loc[
        :,
        [
            "sample",
            "patient",
            "tissue",
            "origin",
            "condition",
            "tumor_stage",
            "dataset",
        ],
    ]
    .assign(sample_lc=lambda x: x["sample"].str.lower())
    .drop_duplicates()
    .set_index("sample_lc")
)

# %%
ad_cpdb.obs = ad_cpdb.obs.join(sample_info)

# %%
ad_cpdb = ad_cpdb[
    :,
    (ad_cpdb.var["cluster_1"] != ad_cpdb.var["cluster_2"])
    & (ad_cpdb.var["cluster_1"] != "other")
    & (ad_cpdb.var["cluster_2"] != "other"),
].copy()

# %%
ad_cpdb.shape

# %% [markdown]
# ## Store results

# %%
ad_cpdb.write_h5ad(f"{artifact_dir}/adata_cpdb.h5ad")

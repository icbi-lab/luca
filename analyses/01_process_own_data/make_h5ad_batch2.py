# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

import scanpy as sc
import pandas as pd
from glob import glob
from pathlib import Path
import re
import scipy.sparse
from multiprocessing import Pool
import anndata

filenames = glob("../../data/11_own_datasets/batch2_5patients/processed/*.csv")

filenames

meta = (
    pd.read_excel(
        "../../tables/patient_table_batch2_5_patients.xlsx", engine="openpyxl"
    )
    .dropna(how="all")
    .dropna(axis="columns", how="all")
)

meta


def load_counts(meta):
    p = Path(f"../../data/11_own_datasets/batch2_5patients/processed/{meta['file_id']}_{'N' if meta['origin'] == 'normal_adjacent' else 'T'}.csv")
    expr = pd.read_csv(p, skiprows=5, index_col=0)
    gene_expr = expr.loc[
        (~expr.index.str.contains("\(Ab\)") & ~expr.index.str.startswith("Lex_")), :
    ]
    obs = pd.DataFrame().assign(cell_id=expr.columns, **meta)
    obs.set_index("cell_id", inplace=True)
    adata = sc.AnnData(X=scipy.sparse.csc_matrix(gene_expr).T, obs=obs)
    adata.var_names = gene_expr.index
    return adata


with Pool(16) as p:
    adatas = p.map(load_counts, meta.to_dict(orient="records"))

adatas[0]

adata = anndata.concat(adatas, index_unique="_", join="outer")

adata.obs["sample"] = [f"{patient}_{origin}" for patient, origin in zip(adata.obs["patient"], adata.obs["origin"])]
adata.obs["platform"] = "BD-Rhapsody"
adata.obs["platform_fine"] = "BD-Rhapsody"

adata.obs.drop_duplicates()

adata

adata.write_h5ad(
    "../../data/11_own_datasets/batch2_5patients/h5ad_raw/ukim_v_batch2.h5ad",
    compression="lzf",
)

# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python [conda env:.conda-scanpy_2020-12]
#     language: python
#     name: conda-env-.conda-scanpy_2020-12-py
# ---

import scanpy as sc
import pandas as pd
from glob import glob
from pathlib import Path
import re
import scipy.sparse
from multiprocessing import Pool
import anndata

filenames = glob("../../data/11_own_datasets/batch1_3patients/processed/*.csv")

meta = pd.read_excel("../../tables/patient_table_batch1_3_patients.xlsx", engine="openpyxl", skiprows=1, skipfooter=1)

meta.rename(columns={"Tumornummer": "tumor_id", "Patient": "patient", "Alter": "age", "Geschlecht": "sex", "Tumor": "tumor_type"}, inplace=True)
meta["sex"] = [{"W": 'f', "M":'m'}[x] for x in meta["sex"]]


def load_counts(filename):
    p = Path(filename)
    ((patient, tissue), ) = re.findall("(P\d+)_(.*)\.csv", p.name)
    expr = pd.read_csv(p, skiprows=5, index_col=0)
    gene_expr = expr.loc[(~expr.index.str.contains("\(Ab\)") & ~expr.index.str.startswith("Lex_")), :]
    ab_expr = expr.loc[expr.index.str.contains("\(Ab\)"), :]
    obs = pd.DataFrame().assign(cell_id=expr.columns)
    obs["patient"] = patient
    obs["tissue"] = tissue
    obs.set_index("cell_id", inplace=True)
    adata = sc.AnnData(X=scipy.sparse.csc_matrix(gene_expr).T, obs=obs)
    adata.var_names = gene_expr.index
    adata.obsm["surface_protein"] = ab_expr.T
    return adata


with Pool(16) as p:
    adatas = p.map(load_counts, filenames)

adata = anndata.concat(adatas, index_unique="_", join="outer")

meta["patient"] = [x.strip() for x in meta["patient"]]

adata.obs = adata.obs.reset_index().merge(meta, on=["patient"], how="left").set_index("cell_id")

adata

adata.obs.drop_duplicates()

adata.obs["condition"] = "NSCLC"
adata.obs["origin"] = ["tumor_primary" if c == "Tumor" else "normal_adjacent" for c in adata.obs["tissue"]]
adata.obs["sample"] = [f"{patient}_{origin}" for patient, origin in zip(adata.obs["patient"], adata.obs["origin"])]
adata.obs["sex"] = [{"m": "male", "f": "female"}[s] for s in adata.obs["sex"]]
adata.obs["tissue"] = "lung"

adata.obs.drop_duplicates()

# !mkdir -p "../../data/11_own_datasets/batch1_3patients/h5ad_raw"

adata.write_h5ad("../../data/11_own_datasets/batch1_3patients/h5ad_raw/batch1_3patients.h5ad", compression="lzf")



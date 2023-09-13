# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python
#     language: python
#     name: python3
# ---

# %%
# %load_ext autoreload
# %autoreload 2


import mnnpy
import anndata
import os

# %%
paths = [file for file in os.listdir(".") if file.endswith(".h5ad")]

if len(paths) != 1:
    raise ValueError("Expected only one file in this directory, but found: {}".format(paths))
adata_path = paths[0]
adata_path

# %%
adata = sc.read_h5ad(adata_path)
adata

# %%
# split by batch for mnnpy
adata_list = []

for batch in set(adata.obs.batch):
    adata_list.append(adata[adata.obs.batch==batch, ].copy())

# %%
# integration with mnn
corrected = mnnpy.mnn_correct(*adata_list, batch_categories= list(set(adata.obs.batch)))

# %%
corrected_values = corrected[0].X

# %%
adata.layers["raw_counts"] = adata.X.copy()
adata.X = corrected_values

# %%
adata.write_h5ad("mnn.h5ad")

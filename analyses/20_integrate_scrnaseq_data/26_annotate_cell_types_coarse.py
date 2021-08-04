# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
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
from scanpy_helpers.annotation import AnnotationHelper
import scvi
import warnings
import numpy as np
from nxfvars import nxfvars

warnings.filterwarnings("ignore", category=FutureWarning)

# %%
ah = AnnotationHelper()

# %%
input_dir = nxfvars.get("input_dir", "../../data/20_integrate_scrnaseq_data/25_merge_solo/")

# %%
adata = sc.read_h5ad(f"{input_dir}/adata.doublet_filtered.h5ad")

# %%

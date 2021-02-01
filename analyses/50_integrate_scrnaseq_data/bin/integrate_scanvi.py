#!/usr/bin/env python
"""Run scANVI integration on prepared anndata object. 

Usage ./integrate_scanpy.py adata.h5ad output_adata.h5ad
"""

import scanpy as sc
import numpy as np
import scvi
import sys

adata_in = sys.argv[1]
adata_out = sys.argv[2]

adata = sc.read_h5ad(adata_in)

sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3", n_top_genes=6000, batch_key="dataset", subset=True
)

adata.obs["batch"] = [
    f"{dataset}_{sample}"
    for dataset, sample in zip(adata.obs["dataset"], adata.obs["sample"])
]

scvi.data.setup_anndata(
    adata,
    batch_key="batch",
    labels_key="cell_type",
)

lvae = scvi.model.SCANVI(adata, "unknown", use_cuda=True, n_latent=20, n_layers=2)

lvae.train(n_epochs_semisupervised=100)

adata.obs["cell_type_predicted"] = lvae.predict(adata)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=10)
sc.tl.umap(adata)

adata.write_h5ad(adata_out)

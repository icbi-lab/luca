#!/usr/bin/env python
"""Run scANVI integration on prepared anndata object.

Usage ./integrate_scanpy.py adata.h5ad output_adata.h5ad model_out_dir
"""

import scanpy as sc
import numpy as np
import scvi
import sys

early_stopping_kwargs = {
    "early_stopping_metric": "elbo",
    "save_best_state_metric": "elbo",
    "patience": 10,
    "threshold": 0,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}
early_stopping_kwargs_scanvi = {
    "early_stopping_metric": "accuracy",
    "save_best_state_metric": "accuracy",
    "on": "full_dataset",
    "patience": 10,
    "threshold": 0.001,
    "reduce_lr_on_plateau": True,
    "lr_patience": 8,
    "lr_factor": 0.1,
}

vae_epochs = 500
scanvi_epochs = 200


adata_in = sys.argv[1]
adata_out = sys.argv[2]
model_out = sys.argv[3]
use_highly_variable = bool(int(sys.argv[4]))

adata = sc.read_h5ad(adata_in)

if use_highly_variable:
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

lvae = scvi.model.SCANVI(adata, "unknown", use_cuda=True, n_latent=30, n_layers=2)

lvae.train(
    n_epochs_unsupervised=vae_epochs,
    n_epochs_semisupervised=scanvi_epochs,
    lr=0.0001,
    unsupervised_trainer_kwargs=dict(early_stopping_kwargs=early_stopping_kwargs),
    semisupervised_trainer_kwargs=dict(
        metrics_to_monitor=["elbo", "accuracy"],
        early_stopping_kwargs=early_stopping_kwargs_scanvi,
    ),
)
lvae.save(model_out)

adata.obs["cell_type_predicted"] = lvae.predict(adata)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=10)
sc.tl.umap(adata)

adata.write_h5ad(adata_out)

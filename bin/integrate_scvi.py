#!/usr/bin/env python

import scanpy as sc
import numpy as np
import scvi
import sys
import argparse


def set_all_seeds(seed=0):
    import os
    import random
    import numpy as np
    import torch

    scvi.settings.seed = seed
    os.environ["PYTHONHASHSEED"] = str(seed)  # Python general
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    np.random.seed(seed)  # Numpy random
    random.seed(seed)  # Python random

    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # For multiGPU


def main(
    adata_in,
    *,
    adata_out,
    model_out,
    use_highly_variable: bool,
    batch_key,
    hvg_batch_key=None,
    labels_key=None,
    use_arches_params=False,
):
    """
    Run scvi tools integration

    Parameters
    ----------
    adata_in
        Input anndata containing all datasets with batch annotation in obs and
        raw counts in X
    adata_out
        Path to write integrated anndata object to (will contain obsm[X_{algorithm}])
    model_out
        Path to write the scVI model to.
    use_highly_variable
        Whether or not to subset to the 6000 most highly variable genes before
        integration
    batch_key
        Key in adata.obs containing the batch annotation
    hvg_batch_key
        Key in adata.obs containing the batch annotation used for selecting
        highly variable genes. May be different from batch_key, as highly_variable_genes
        might have trouble to work with small batches (e.g. choose dataset instead
        of sample).
    labels_key
        Key in adata.obs containing cell-type labels for semi-supervised training
        and later SCANVI inference
    use_arches_params
        Use parameters that were found to work well with scArches, cf. https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scarches_scvi_tools.html#Online-update-of-SCANVI
    """
    set_all_seeds()

    adata = sc.read_h5ad(adata_in)

    extra_params = {}
    if use_arches_params:
        extra_params = dict(
            use_layer_norm="both",
            use_batch_norm="none",
            encode_covariates=True,
            dropout_rate=0.2,
            n_layers=2,
        )

    if use_highly_variable:
        if hvg_batch_key is None:
            hvg_batch_key = batch_key
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=6000,
            batch_key=hvg_batch_key,
            subset=True,
        )

    def _setup_anndata(*args, **kwargs):
        """Setup anndata compatible with different scvi-tools versions"""
        try:
            return scvi.data.setup_anndata(*args, **kwargs)
        except AttributeError:
            return scvi.model.SCVI.setup_anndata(*args, **kwargs)

    if labels_key is not None:
        _setup_anndata(adata, batch_key=batch_key, labels_key=labels_key)
    else:
        _setup_anndata(adata, batch_key=batch_key)

    vae = scvi.model.SCVI(adata, **extra_params)

    vae.train(
        # TODO: check how we can use the GPU
        # use_gpu=True,
    )
    vae.save(model_out)

    adata.obsm[f"X_scVI"] = vae.get_latent_representation(adata)

    adata.write_h5ad(adata_out)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Integrate scRNA-seq data with scVI")
    parser.add_argument(
        "adata_in",
        type=str,
        help="Anndata containing all samples/datasets to integrate",
    )
    parser.add_argument(
        "--adata_out",
        dest="adata_out",
        type=str,
        help="Output path for integrated anndata object",
    )
    parser.add_argument(
        "--model_out",
        dest="model_out",
        type=str,
        help="Output path for the trained scVI model (directory)",
    )
    parser.add_argument(
        "--use_hvg",
        dest="use_hvg",
        type=bool,
        help="Whether to use all genes (False) or the 6000 most highly variable genes (True)",
        default="1",
    )
    parser.add_argument(
        "--batch_key",
        dest="batch_key",
        type=str,
        help="Key in adata.obs that contains the batch annotation",
        default="batch",
    )
    parser.add_argument(
        "--hvg_batch_key",
        dest="hvg_batch_key",
        type=str,
        help="Key in adata.obs that contains the batch annotation used for highly variable genes",
        default="batch",
    )
    parser.add_argument(
        "--labels_key",
        dest="labels_key",
        type=str,
        help="Key in adata.obs that contains the cell-type labels (optional)",
        default=None,
    )
    parser.add_argument(
        "--arches_params",
        dest="arches_params",
        type=bool,
        help="Run with parameters optimized for scArches mode",
        default=False,
    )
    args = parser.parse_args()
    main(
        args.adata_in,
        adata_out=args.adata_out,
        model_out=args.model_out,
        use_highly_variable=args.use_hvg,
        batch_key=args.batch_key,
        hvg_batch_key=args.hvg_batch_key,
        labels_key=args.labels_key,
        use_arches_params=args.arches_params,
    )

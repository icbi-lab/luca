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
    algorithm="scVI",
):
    set_all_seeds()

    adata = sc.read_h5ad(adata_in)

    if use_highly_variable:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=6000,
            batch_key=batch_key,
            subset=True,
        )

    scvi.data.setup_anndata(
        adata,
        batch_key=batch_key,
    )

    # lvae = scvi.model.SCANVI(adata, "unknown", use_cuda=True, n_latent=30, n_layers=2)
    lvae = scvi.model.SCVI(adata)

    lvae.train(
        use_gpu=True,
    )
    lvae.save(model_out)

    # adata.obs["cell_type_predicted"] = lvae.predict(adata)
    adata.obsm[f"X_{algorithm}"] = lvae.get_latent_representation(adata)

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
        help="Whether to use all genes (False) or the 6000 most highly variable genese (True)",
        default="1",
    )
    parser.add_argument(
        "--batch_key",
        dest="batch_key",
        type=str,
        help="Key in adata.obs that contains the batch annotation",
        default="batch",
    )
    args = parser.parse_args()
    main(
        args.adata_in,
        adata_out=args.adata_out,
        model_out=args.model_out,
        use_highly_variable=args.use_hvg,
        batch_key=args.batch_key,
    )

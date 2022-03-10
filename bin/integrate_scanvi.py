#!/usr/bin/env python

import scanpy as sc
import numpy as np
import scvi
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
    model_in,
    *,
    adata_out,
    model_out,
    batch_key="batch",
    labels_key="cell_type",
):
    """
    Run scANVI on a pre-trained SCVI model

    Parameters
    ----------
    adata_in
        Input anndata containing all datasets with batch annotation in obs and
        raw counts in X
    model_in
        pre-trained scVI model
    adata_out
        ouput anndata (contains adata.obsm["scANVI"])
    model_out
        Output directory for the scANVI model
    labels_key
        Key in adata.obs which contains the cell-type labels. Unlabelled cells
        must have the label "unknown"
    """
    set_all_seeds()

    print(f"batch_key={batch_key}, labels_key={labels_key}")

    def _setup_lvae(adata, vae, batch_key, labels_key):
        """Setup anndata compatible with different scvi-tools versions"""
        try:
            scvi.data.setup_anndata(adata, batch_key=batch_key, labels_key=labels_key)
            lvae = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="unknown", adata=adata)
        except AttributeError:
            scvi.model.SCANVI.setup_anndata(adata, unlabeled_category="unknown", batch_key=batch_key, labels_key=labels_key)
            lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=adata, unlabeled_category="unknown", labels_key=labels_key)
        return lvae

    adata = sc.read_h5ad(adata_in)

    vae = scvi.model.SCVI.load(model_in, adata, use_gpu=False)
    lvae = _setup_lvae(adata, vae, batch_key=batch_key, labels_key=labels_key)

    lvae.train(use_gpu=False)
    lvae.save(model_out)

    adata.obsm[f"X_scANVI"] = lvae.get_latent_representation(adata)
    adata.obs[f"{labels_key}_predicted"] = lvae.predict(adata)

    adata.write_h5ad(adata_out)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Integrate scRNA-seq data with scANVI")
    parser.add_argument(
        "adata_in",
        type=str,
        help="Anndata containing all samples/datasets to integrate",
    )
    parser.add_argument(
        "model_in",
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
        "--batch_key",
        dest="batch_key",
        type=str,
        help="Key in adata.obs that contains the batch annotation",
        default="batch",
    )
    parser.add_argument(
        "--labels_key",
        dest="labels_key",
        type=str,
        help="Key in adata.obs that contains the cell-type annotation",
        default="cell_type",
    )
    args = parser.parse_args()
    main(
        args.adata_in,
        args.model_in,
        adata_out=args.adata_out,
        model_out=args.model_out,
        batch_key=args.batch_key,
        labels_key=args.labels_key,
    )

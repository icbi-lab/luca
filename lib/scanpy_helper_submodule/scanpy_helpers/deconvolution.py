import numpy as np
from anndata import AnnData


def balanced_subsample(
    adata: AnnData,
    *,
    cell_type_key: str,
    patient_key: str,
    n_each: int = 10,
    seed: int = 0
) -> AnnData:
    """Subsample a dataset such that cell-types and patients are balanced.

    This helps generating a signature matrix that is not biased towards a specific patient.

    Subsampling *without replacement*, i.e. if less than `n_each` cells are in a group,
    only those cells will be returned and not over-sampled.

    Parameters
    ----------
    adata
        anndata object with cells of multiple cell-types from multiple patients
    cell_type_key
        column in adata.obs that holds the cell-type annotation
    patient_key
        columnin adata.obs that holds the patient (or batch) information
    n_each
        Maximum number of cells to sample from each patient for each cell_type
    seed
        random seed for reproducibility

    Returns
    -------
    Subset on `adata`
    """
    np.random.seed(seed)
    patients = adata.obs[patient_key].unique()
    cell_types = adata.obs[cell_type_key].unique()
    cell_idx = []
    for patient in patients:
        mask_patient = adata.obs[patient_key] == patient
        for cell_type in cell_types:
            mask = mask_patient & (adata.obs[cell_type_key] == cell_type)
            current_idx = adata.obs_names[mask]
            cell_idx.extend(
                np.random.choice(
                    current_idx, size=min(n_each, np.sum(mask)), replace=False
                )
            )

    return adata[cell_idx, :]

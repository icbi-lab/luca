import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse
import mygene
import itertools
import warnings

"""Definitions for metadata consistency checks"""
MANDATORY_COLS = [
    "sample",
    "patient",
    "tissue",
    "origin",
    "condition",
    "dataset",
    "sex",
]

VALID_ORIGIN = [
    "tumor_primary",
    "normal_adjacent",
    "tumor_edge",
    "tumor_metastasis",
    "blood_peripheral",
    "lymph_node",
    "normal",
    "nan",
]

VALID_TISSUE = ["lung", "adrenal", "brain", "liver", "lymph_node", "pleura"]

VALID_CONDITION = ["LUAD", "NSCLC", "LSCC", "COPD", "healthy_control"]

VALID_SEX = ["male", "female", "nan"]


def sanitize_adata(adata):
    # sample should be dataset specific
    adata.obs["sample"] = [
        f"{dataset}_{sample}"
        for dataset, sample in zip(adata.obs["dataset"], adata.obs["sample"])
    ]

    # X should be in CSR format
    adata.X = adata.X.tocsr()
    adata._sanitize()


def drop_duplicated_genes(adata):
    duplicated = adata.var_names.duplicated()
    if np.sum(duplicated):
        warnings.warn(f"Removing {np.sum(duplicated)} genes. ")
        return adata[:, ~duplicated].copy()
    else:
        return adata


def validate_adata(adata):
    _validate_obs(adata)

    assert isinstance(adata.X, scipy.sparse.csr_matrix)

    assert adata.obs_names.is_unique, "Obs names not unique"
    assert adata.var_names.is_unique, "var names not unique"

    # X should be all integers
    assert np.all(np.modf(adata.X.data)[0] == 0), "X does not contain all integers"


def _validate_obs(adata):
    """Consistency checks for adata.obs"""

    obs = adata.obs

    for col in MANDATORY_COLS:
        assert col in obs.columns, "{} is a mandatory column".format(col)

    # Only one patient per sample
    sample_count = obs.groupby(["sample", "patient"], observed=True)["sample"].nunique()
    assert np.all(
        sample_count == 1
    ), "sample must be unique for each patient, origin and replicate"

    def _check_col(col, valid_keywords):
        isin_col = obs[col].isin(valid_keywords)
        assert np.all(
            isin_col
        ), f"Invalid words in {col}: {np.unique(obs[col].values[~isin_col])}"

    # check controlled vocabulary
    _check_col("origin", VALID_ORIGIN)
    _check_col("tissue", VALID_TISSUE)
    _check_col("condition", VALID_CONDITION)
    _check_col("sex", VALID_SEX)


def undo_log_norm(adata):
    """Reverse a log-normalization, assuming that each
    sample has at least one cell with exactely one count. This assumption
    is reasonable, at least for 10x data. (I checked on the
    Lambrechts dataset and it is true there).
    """
    adata.layers["log_norm"] = adata.X.copy()
    x_log_norm = adata.X.tocsr()
    x_log_norm.sort_indices()
    x_norm = x_log_norm.copy()
    x_norm.data = np.expm1(x_norm.data)

    # assuming that each sample has at least one cell with exactely one count
    size_factors = np.array([np.min(x_norm[i, :].data) for i in range(x_norm.shape[0])])
    x_raw_counts = scipy.sparse.diags(1 / size_factors) @ x_norm
    x_raw_counts.data = np.rint(x_raw_counts.data)
    x_raw_counts.sort_indices()
    adata.X = x_raw_counts.copy()


def normalize_by_gene_length(adata) -> sc.AnnData:
    """Normalize Smart-seq2 data by gene length for scANVI.

    This follows the code in the scANVI tutorial.

    Assumes that there is a layer "counts_length_scaled" with lengh-normalized
    gene counts and a that raw counts are in X.
    """
    gene_lengths = adata.layers["counts_length_scaled"].data / adata.X.data
    median_gene_length = np.median(1 / gene_lengths)
    x_length_scaled = adata.layers["counts_length_scaled"] * median_gene_length
    x_length_scaled.data = np.rint(x_length_scaled.data)
    return sc.AnnData(
        X=x_length_scaled,
        var=adata.var,
        obs=adata.obs,
    )


def add_doublet_annotation(adata, doublet_file, plot_title):
    """Add the doublet annotation based on the doublet file generated
    by the single-cell-analysis-nf pipeline.

    Returns a temporary anndata object with umap visualized"""
    doublets = pd.read_csv(
        doublet_file,
        names=["cell_id", "is_doublet"],
    ).set_index("cell_id")
    adata.obs["is_doublet"] = doublets.reindex(adata.obs_names)["is_doublet"].astype(
        "str"
    )
    adata_vis = adata.copy()
    sc.pp.normalize_total(adata_vis)
    sc.pp.log1p(adata_vis)

    sc.pp.highly_variable_genes(adata_vis, n_top_genes=4000, flavor="cell_ranger")
    sc.tl.pca(adata_vis)
    sc.pp.neighbors(adata_vis)
    sc.tl.umap(adata_vis)
    sc.pl.umap(adata_vis, color=["sample", "is_doublet"], title=plot_title)
    return adata_vis


def remap_gene_symbols(adata):
    mg = mygene.MyGeneInfo()
    query_result = mg.querymany(
        adata.var_names.values,
        scopes="symbol,name,alias",
        fields="symbol",
        species="human",
    )
    query_result = sorted(
        query_result, key=lambda r: (r["query"], r.get("_score", float("inf")))
    )

    # Take the result with the highest score when multiple hits are found.
    # Take the query if no matching symbol was found
    gene_dict = {
        query: next(results).get("symbol", query)
        for query, results in itertools.groupby(query_result, key=lambda r: r["query"])
    }

    adata.var["original_gene_symbol"] = adata.var_names
    adata.var_names = [gene_dict[g] for g in adata.var_names.values]

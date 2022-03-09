from scanpy_helpers.integration import aggregate_duplicate_gene_symbols
import scanpy as sc
import numpy as np
import pandas as pd
import numpy.testing as npt
import scipy.sparse


def test_aggregate_duplicate_gene_symbols():
    adata = sc.AnnData(
        X=scipy.sparse.csc_matrix([[4, 5, 1], [1, 6, 1]]),
        var=pd.DataFrame().assign(idx=["s", "s", "x"]).set_index("idx"),
    )
    adata_dedup = aggregate_duplicate_gene_symbols(adata)
    npt.assert_equal(adata_dedup.X.toarray(), np.array([[5, 1], [6, 1]]))
    npt.assert_equal(adata_dedup.var_names.values, ["s", "x"])

from scanpy_helpers.integration import aggregate_duplicate_gene_symbols
import scanpy as sc
import numpy as np
import pandas as pd
import numpy.testing as npt


def test_aggregate_duplicate_gene_symbols():
    adata = sc.AnnData(
        X=np.array([[4, 5, 1], [1, 6, 1]]),
        var=pd.DataFrame().assign(idx=["s", "s", "x"]).set_index("idx"),
    )
    adata_dedup = aggregate_duplicate_gene_symbols(adata)
    npt.assert_equal(adata_dedup.X, np.array([[5, 1], [6, 1]]))
    npt.assert_equal(adata_dedup.var_names.values, ["s", "x"])

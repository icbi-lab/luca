from anndata import AnnData
import numpy as np
import pandas as pd
import scipy.sparse as sp
import pytest


@pytest.fixture(params=[np.array, sp.csr_matrix, sp.csc_matrix])
def adata_hierarchical1(request):
    return AnnData(
        X=request.param(
            np.array(
                [
                    [10, 12, 1, 10, 3, 20, 13, 5],
                    [0, 0, 0, 0, 0, 1, 0, 0],
                ]
            ).T
        ),
        obs=pd.DataFrame(index=["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"]).assign(
            dataset=["d1", "d1", "d2", "d1", "d2", "d3", "d1", "d2"],
            patient=["p1", "p2", "p1", "p1", "p1", "p1", "p3", "p2"],
        ),
        var=pd.DataFrame(index=["A", "B"]),
    )


@pytest.fixture(params=[np.array, sp.csr_matrix, sp.csc_matrix])
def adata_ithgex(request):
    return AnnData(
        X=request.param(
            np.array(
                [
                    [1, 1, 1, 1, 1, 1, 2, 3],
                    [2, 2, 2, 2, 2, 2, 8, 0],
                    [3, 3, 3, 3, 3, 10, 3, 7],
                ]
            ).T
        ),
        obsm={
            "X_cnv": request.param(
                np.array(
                    [
                        [1, 1, 1, 2, 2, 1, 1, 1],
                        [2, 2, 2, 1, 1, 2, 2, 2],
                        [4, 4, 4, 2, 2, 3, 3, 3],
                        [2, 2, 2, 4, 4, 4, 4, 4],
                    ]
                ).T
            )
        },
        obs=pd.DataFrame(index=["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"]).assign(
            group=list("AAAAABBB"),
        ),
        var=pd.DataFrame(index=["x", "y", "z"]),
    )


# @pytest.fixture(params=[np.array, sp.csr_matrix, sp.csc_matrix])
# def adata_multi_group(request):
#     """Adata with multiple groups which don't follow a hierarchical structure."""
#     return AnnData(
#         X=request.param(
#             np.array(
#                 [
#                     [10, 12, 1, 10, 3, 20, 13, 5],
#                     [0, 0, 0, 0, 0, 1, 0, 0],
#                 ]
#             ).T
#         ),
#         obs=pd.DataFrame(index=["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"]).assign(
#             patient=["p1", "p2", "p1", "p1", "p1", "p1", "p3", "p2"],
#             origin=["T", "T", "T", "N", "N", "T", "T", "N"],
#         ),
#         var=pd.DataFrame(index=["A", "B"]),
#     )

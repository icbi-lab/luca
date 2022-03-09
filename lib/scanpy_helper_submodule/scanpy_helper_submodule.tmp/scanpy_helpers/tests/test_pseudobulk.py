from .fixtures import adata_hierarchical1
from scanpy_helpers.pseudobulk import pseudobulk
import pytest
import numpy.testing as npt
import numpy as np
import pandas as pd
import pandas.testing as pdt


@pytest.mark.parametrize(
    "groupby,min_obs,obs_expected,x_expected,aggr_fun",
    [
        # np.mean tests
        (
            ["dataset", "patient"],
            0,
            [
                ["d1", "p1", 2],
                ["d1", "p2", 1],
                ["d2", "p1", 2],
                ["d3", "p1", 1],
                ["d1", "p3", 1],
                ["d2", "p2", 1],
            ],
            [[10, 0], [12, 0], [2, 0], [20, 1], [13, 0], [5, 0]],
            np.mean,
        ),
        (
            ["dataset", "patient"],
            2,
            [
                ["d1", "p1", 2],
                ["d2", "p1", 2],
            ],
            [[10, 0], [2, 0]],
            np.mean,
        ),
        (
            "dataset",
            1,
            [
                ["d1", 4],
                ["d2", 3],
                ["d3", 1],
            ],
            [[11.25, 0], [3, 0], [20, 1]],
            np.mean,
        ),
        # np.sum tests
        (
            ["dataset", "patient"],
            0,
            [
                ["d1", "p1", 2],
                ["d1", "p2", 1],
                ["d2", "p1", 2],
                ["d3", "p1", 1],
                ["d1", "p3", 1],
                ["d2", "p2", 1],
            ],
            [[20, 0], [12, 0], [4, 0], [20, 1], [13, 0], [5, 0]],
            np.sum,
        ),
        (
            ["dataset", "patient"],
            2,
            [
                ["d1", "p1", 2],
                ["d2", "p1", 2],
            ],
            [[20, 0], [4, 0]],
            np.sum,
        ),
        (
            "dataset",
            1,
            [
                ["d1", 4],
                ["d2", 3],
                ["d3", 1],
            ],
            [[45, 0], [9, 0], [20, 1]],
            np.sum,
        ),
    ],
)
def test_pseudobulk1(
    adata_hierarchical1, groupby, obs_expected, x_expected, min_obs, aggr_fun
):
    x_expected = np.array(x_expected)
    obs_expected = pd.DataFrame(
        obs_expected,
        columns=([groupby] if isinstance(groupby, str) else groupby) + ["n_obs"],
    )
    obs_expected.index = obs_expected.index.astype(str)
    adata_pb = pseudobulk(
        adata_hierarchical1, groupby=groupby, min_obs=min_obs, aggr_fun=aggr_fun
    )
    pdt.assert_frame_equal(adata_pb.obs, obs_expected)
    npt.assert_almost_equal(adata_pb.X, x_expected)

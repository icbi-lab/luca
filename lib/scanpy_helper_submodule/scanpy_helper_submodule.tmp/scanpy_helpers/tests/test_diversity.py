from .fixtures import adata_ithgex
import pytest
from scanpy_helpers.diversity import ithgex, ithcna


def test_ithgex(adata_ithgex):
    res = ithgex(adata_ithgex, "group", inplace=False)
    assert res["A"] == 0
    assert res["B"] == pytest.approx(1.2628, abs=0.001)


def test_ithcna(adata_ithgex):
    res = ithcna(adata_ithgex, "group", inplace=False)
    assert res["A"] == pytest.approx(1.053, abs=0.001)
    assert res["B"] == 0

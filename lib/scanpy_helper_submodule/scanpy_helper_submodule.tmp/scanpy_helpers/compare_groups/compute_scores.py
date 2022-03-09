"""Run dorothea, progeny and cytosig with parameters suitable for between-sample comparisons.

In particular, permutation tests are not necessary, as we perform statistics on the sample-level.
"""

from functools import lru_cache
from ..util import suppress_stdout
import pandas as pd
from .. import assets
import importlib.resources as pkg_resources


@lru_cache
def _progeny_model():
    import progeny

    return progeny.load_model(
        organism="Human",  # If working with mouse, set to Mouse
        top=1000,  # For sc we recommend ~1k target genes since there are dropouts
    )


@lru_cache
def _dorothea_model():
    import dorothea

    return dorothea.load_regulons(
        [
            "A",
            "B",
        ],  # Which levels of confidence to use (A most confident, E least confident)
        organism="Human",  # If working with mouse, set to Mouse
    )


@lru_cache
def _cytosig_model():
    return pd.read_csv(
        pkg_resources.open_text(assets, "cytosig_signature_matrix.tsv"),
        sep="\t",
        index_col=0,
    )


@suppress_stdout
def run_progeny(adata):
    import progeny

    tmp_adata = adata.copy()

    progeny.run(
        tmp_adata,  # Data to use
        _progeny_model(),  # PROGENy network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # Pathways with less than 5 targets will be ignored
    )
    return progeny.extract(tmp_adata)


@suppress_stdout
def run_dorothea(adata):
    import dorothea

    tmp_adata = adata.copy()
    dorothea.run(
        tmp_adata,  # Data to use
        _dorothea_model(),  # Dorothea network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # TF with less than 5 targets will be ignored
    )
    return dorothea.extract(tmp_adata)


@suppress_stdout
def run_cytosig(adata):
    """Run cytosig.

    The algorithm in dorothea-py gives the same results as the original implementation (it's a simple matrix multiplication)
    and is a lot faster.
    """
    import progeny

    tmp_adata = adata.copy()
    progeny.run(
        tmp_adata,  # Data to use
        _cytosig_model(),  # PROGENy network
        center=True,  # Center gene expression by mean per cell
        num_perm=0,  # Simulate m random activities
        norm=True,  # Normalize by number of edges to correct for large regulons
        scale=True,  # Scale values per feature so that values can be compared across cells
        use_raw=True,  # Use raw adata, where we have the lognorm gene expression
        min_size=5,  # Pathways with less than 5 targets will be ignored
        obsm_key="cytosig",
    )
    return progeny.extract(tmp_adata, obsm_key="cytosig")

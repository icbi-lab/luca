"""High-level wrappers around linear models to find differences between groups. """


from typing import Dict, List, Mapping, Optional, Sequence, Tuple
from unittest.util import strclass
import warnings
import pandas as pd
from patsy import PatsyError
from tqdm.auto import tqdm
import statsmodels.formula.api as smf
from tqdm.contrib.concurrent import process_map
import re
from anndata import AnnData
from ..pseudobulk import pseudobulk
import numpy as np
import pandas as pd
from ..util import chunk_adatas, fdr_correction
import itertools
from multiprocessing import cpu_count
from threadpoolctl import threadpool_limits
from statsmodels.regression.linear_model import RegressionResultsWrapper


def lm_test_all(
    adatas: Mapping[str, AnnData],
    *,
    groupby: List[str],
    column_to_test: str,
    contrasts: str = "Sum",
    lm_covariate_str: str = "",
    min_categories: int = 2,
    **kwargs,
) -> Optional[pd.DataFrame]:
    """High-level function to compare groups multiple single cell objects.

    Performs pseudobulk aggregation and a statistical test with linear model

    Parameters
    ----------
    adatas
        Dictionary id -> adata
    groupby
        Grouping variables, will be used to generate pseudobulk (e.g. ["patient", "dataset"])
    column_to_test
        column in obs that contains the variables to test. Must be the same within each group.
        (e.g. a property that applies to all cells of a patient)
    contrasts
        Which contrast-coding to use. T = treatment-coding, S =sum-to-zero
        See https://www.statsmodels.org/devel/contrasts.html
    lm_covariate_str
        Covariate string that is literally included in the linear model formula.
        E.g. `+ dataset`.
    min_categories
        Only perform test if there are at least `min_categories` unique entries in `column_to_test`.

    Returns
    -------
    Data frame with coefficients and pvalues.
    """
    res_list = []
    for ct, tmp_adata in tqdm(adatas.items()):
        tmp_bdata = pseudobulk(
            tmp_adata,
            groupby=groupby + [column_to_test],
            aggr_fun=np.mean,
        )
        if tmp_bdata.obs[column_to_test].nunique() >= min_categories:
            try:
                tmp_res = test_lm(
                    tmp_bdata,
                    f"~ C({column_to_test}, {contrasts}) {lm_covariate_str}",
                    column_to_test,
                    contrasts=contrasts,
                    progress=False,
                    **kwargs,
                )
                res_list.append(tmp_res.assign(cell_type=ct))
            except PatsyError as e:
                warnings.warn(f"Iteration on {ct} failed with {e}. Ignoring error. ")

    return pd.concat(res_list).pipe(fdr_correction).sort_values("pvalue")


def test_lm(
    pseudobulk: AnnData,
    formula: str,
    groupby: str,
    *,
    contrasts: str = "Sum",
    robust: bool = False,
    progress: bool = True,
    n_jobs: int = None,
    chunksize=200,
):
    """
    Use a linear model to find differences between groups

    Parameters
    ----------
    pseudobulk
        AnnData object with pseudobulk
    formula
        Formula for the linear model. Currently MUST have the following structure
        `~ 0 + C({column}, {contrasts}) + ...`.
    groupby
        Column in adata that contains the groups to test
    contrasts
        a patsy contrast included in to linear model formula. May be `Sum`, `Treatment`, or `Treatment('<base level>')`.
        See https://www.statsmodels.org/devel/contrasts.html for more details.
    robust
        If true, use heteroskedasticity adjustment with HC3.
        HC3 has been shown to be superior over the default HC1 in http://datacolada.org/99.
        See also: https://www.statsmodels.org/devel/generated/statsmodels.regression.linear_model.OLSResults.get_robustcov_results.html
    progress
        Show the tqdm progress bar
    n_jobs
        Run test in parallel. Set to 1 to disable parallelism.
    chunksize
        Splits up anndata in chunks of var and processes chunks in parallel.

    Returns
    -------
    Pandas data frame with coefficients and pvalues

    """
    if n_jobs is None:
        n_jobs = cpu_count()

    # likely overhead is larger until several times chunksize. Exact optimum not tested.
    if pseudobulk.shape[1] < chunksize * 2 or n_jobs == 1:
        return _test_lm(
            pseudobulk,
            formula,
            groupby,
            contrasts=contrasts,
            robust=robust,
            progress=progress,
        )
    else:
        return pd.concat(
            process_map(
                _test_lm,
                list(chunk_adatas(pseudobulk, chunksize=chunksize)),
                itertools.repeat(formula),
                itertools.repeat(groupby),
                itertools.repeat(contrasts),
                itertools.repeat(False),
                itertools.repeat(robust),
                max_workers=n_jobs,
            )
        )


def _test_all_params(
    res, *, groupby: str, all_groups: Sequence[str], contrasts: str
) -> Tuple[Dict, Dict, Dict]:
    """
    Calculate pvalues and coefficients from a statsmodels results object.

    If we use sum-to-zero coding, perform the test for the omitted level.
    (By default, one level is omitted, because it is redundant. We can compute
    the corresponding pvalue and coefficient based on the others. )

    Parameters
    ----------
    res
        a statsmodels result object coming from a linear model
    groupby
        Name of the column to test
    all_groups
        List of unique values in the `groupby` column
    contrasts
        a patsy contrast included in to linear model formula. May be `Sum`, `Treatment`, or `Treatment('<base level>')`.
        See https://www.statsmodels.org/devel/contrasts.html for more details.

    Returns
    -------
    coefs
        Dictionary mapping keys from the linear model summary data frame onto linear model coefficients
    pvals
        Dictionary mapping keys from the linear model summary data frame onto pvalues (f-test)
    intercept
        Dictionary mapping keys from the linear model summary data frame onto linear model intercepts
    """
    if contrasts.startswith("Sum"):
        contrasts_mode = "sum-to-zero"
    elif contrasts.startswith("Treatment"):
        contrasts_mode = "treatment-coding"
    else:
        raise ValueError(f"Unsupported contrast type: {contrasts}!")

    # in case of treatment coding, we do not want to include the baseline level.
    groups_to_test = (
        all_groups
        if contrasts_mode == "sum-to-zero"
        # the contrast should look like `Treatment("something")`
        else [g for g in all_groups if f"{g}" not in contrasts.replace("Treatment", "")]
    )
    assert "nan" not in groups_to_test

    # Generate keys as they are used in the linear model results data frame
    keys = [f"C({groupby}, {contrasts})[{contrasts[0]}.{g}]" for g in groups_to_test]

    intercept = res.params["Intercept"]

    if contrasts_mode == "sum-to-zero":
        coefs = res.params[keys[:-1]].to_dict()
        pvals = res.pvalues[keys[:-1]].to_dict()

        # test the level that was omitted for redundancy
        coefs[keys[-1]] = -sum(coefs.values())
        pvals[keys[-1]] = float(
            res.f_test(" + ".join([f"{k}" for k in keys[:-1]]) + " = 0").pvalue
        )
    else:
        coefs = res.params[keys].to_dict()
        pvals = res.pvalues[keys].to_dict()

    return coefs, pvals, intercept


def _test_lm(
    pseudobulk: AnnData,
    formula: str,
    groupby: str,
    contrasts: str,
    progress: bool = False,
    robust: bool = False,
) -> pd.DataFrame:
    """Internal helper function for :func:`test_lm` that gets called in parallel."""
    var_names = pseudobulk.var_names

    # Q() contains the dependent variable which will get inserted in each iteration with a `.format()` call.
    formula = "Q('{col}') " + formula

    # add `adata.X` to the dataframe (used by the linear model)
    df = pseudobulk.obs.join(
        pd.DataFrame(pseudobulk.X, index=pseudobulk.obs_names, columns=var_names)  # type: ignore
    )

    # convert to categorical to make sure the levels have an inherent order
    pseudobulk.obs[groupby] = pd.Categorical(pseudobulk.obs[groupby])  # type: ignore
    # only using the categories (rather than just .unique() gets rid of nans)
    all_groups = pseudobulk.obs[groupby].cat.categories.tolist()

    results: List[pd.DataFrame] = []
    var_iter = tqdm(var_names) if progress else var_names
    for col in var_iter:
        try:
            with threadpool_limits(1):
                mod = smf.ols(formula=formula.format(col=col), data=df)
                res = mod.fit()
                if robust:
                    res = RegressionResultsWrapper(
                        res.get_robustcov_results(cov_type="HC3")
                    )
                coefs, pvals, intercept = _test_all_params(
                    res, groupby=groupby, all_groups=all_groups, contrasts=contrasts
                )
                res_df = (
                    pd.DataFrame.from_dict(coefs, orient="index", columns=["coef"])
                    .assign(intercept=intercept)
                    .join(
                        pd.DataFrame.from_dict(
                            pvals, orient="index", columns=["pvalue"]
                        )
                    )
                    .assign(
                        variable=col,
                        # extract plain group variable from model summary data frame:
                        group=lambda x: [
                            re.search(f"\[{contrasts[0]}\.(.*)\]", k).groups()[0]
                            for k in x.index
                        ],
                    )
                )
            results.append(res_df)
        except ValueError:
            pass

    try:
        return pd.concat(results)
    except ValueError:
        return pd.DataFrame()

"""
Functions to find and test gene signatures.

Inspired by MCPCounter and BioQC.

The original cutoffs used in the MCP counter paper are
 * fold change >= ``min_fc`` (default=2)
 * specific fold change >= ``min_sfc`` (default=1.5)
 * AUC ROC >= ``min_auc`` (default=0.97)
"""

from typing import List
import numpy as np
from sklearn.metrics import roc_auc_score
import sklearn.model_selection
from numba import njit
from anndata import AnnData
import itertools
import scipy.stats
from .pseudobulk import pseudobulk
import scanpy as sc
from tqdm.auto import tqdm


def fold_change(
    adata, *, obs_col, positive_class, key_added="fold_change", inplace=True
):
    """
    Compute the fold change between the
    positive and negative samples for a single gene `G`.

    According to `Becht et al.`_ the fold change is defined as
    .. math::
        FC = X - \overline{X}
    where :math:`X` is the mean of positive and :math:`\overline{X}` is the mean of the negative samples.

    Parameters
    ----------
    adata
        annotated data matrix. Values in X are expected to be log-transformed and normalized.
    obs_col
        Column in adata.obs with the annotation to compare
    positive_class
        Value in `adata.obs[obs_col]` to be considered the positive class
    key_added
        Key under which the result will be stored in `adata.var`
    inplace
        If True, store the result in adata.var. Otherwise return result.

    Returns
    -------
    fold change

    .. _Becht et al.:
        http://dx.doi.org/10.1186/s13059-016-1070-5
    """
    positive_mask = adata.obs[obs_col] == positive_class
    fold_change = np.mean(adata.X[positive_mask, :], axis=0) - np.mean(
        adata.X[~positive_mask, :], axis=0
    )

    if inplace:
        adata.var[key_added] = fold_change
    else:
        return fold_change


def specific_fold_change(
    adata, *, obs_col, positive_class, key_added="specific_fold_change", inplace=True
):
    """
    Compute the specific fold change of the positive class with respect to all other classes
    for a single gene `G`.
    According to `Becht et al.`_ the specific fold change is defined as
    .. math::
        sFC = (X - \overline{X}_{min})/(\overline{X}_{max} - \overline{X}_{min})
    where :math:`X` is the mean of positive and :math:`\overline{X}_{max}` is the maximum mean over
    all negative classes and :math:`\overline{X}_{min}` is the minimal mean over all negative classes.

    Parameters
    ----------
    adata
        annotated data matrix. Values in X are expected to be log-transformed and normalized
    obs_col
        Column in adata.obs with the annotation to compare
    positive_class
        Value in `adata.obs[obs_col]` to be considered the positive class
    key_added
        Key under which the result will be stored in `adata.var`
    inplace
        If True, store the result in adata.var. Otherwise return result.

    Returns
    -------
    specific fold change

    .. _Becht et al.:
        http://dx.doi.org/10.1186/s13059-016-1070-5
    """
    neg_classes = [x for x in adata.obs[obs_col].unique() if x != positive_class]
    mean_per_class = np.hstack(
        [
            np.mean(adata.X[adata.obs[obs_col] == tmp_class, :], axis=0)[:, np.newaxis]
            for tmp_class in neg_classes
        ]
    )
    x_min = np.min(mean_per_class, axis=1)
    x_max = np.max(mean_per_class, axis=1)
    spec_fc = np.divide(
        np.mean(adata.X[adata.obs[obs_col] == positive_class, :], axis=0) - x_min,
        x_max - x_min,
    )

    if inplace:
        adata.var[key_added] = spec_fc
    else:
        return spec_fc


def roc_auc(
    adata,
    *,
    obs_col,
    positive_class,
    key_added="roc_auc",
    inplace=True,
    multi_class="ovo"
):
    """
    Compute the Receiver Operator Characteristics Area under the Curve (ROC AUC) for a single gene `G`.
    This tells how well the gene discriminates between the two classes.
    This is a wrapper for the scikit-learn `roc_auc_score`_ function.

    Parameters
    ----------
    adata
        annotated data matrix. Values in X are expected to be log-transformed and normalized
    obs_col
        Column in adata.obs with the annotation to compare
    positive_class
        Value in `adata.obs[obs_col]` to be considered the positive class
    key_added
        Key under which the result will be stored in `adata.var`
    inplace
        If True, store the result in adata.var. Otherwise return result.

    Returns
    -------
    roc auc score

    .. _roc_auc_score:
        http://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score
    >>> expr = np.array([2, 3, 5, 4, 9, 15])
    >>> target = np.array(["A", "B", "C", "A", "B", "C"])
    >>> roc_auc(expr, target == "A")
    0.125
    """

    def _roc_auc(arr, mask):
        return roc_auc_score(mask, arr)

    positive_mask = adata.obs[obs_col] == positive_class
    roc_auc = np.apply_along_axis(_roc_auc, 0, adata.X, positive_mask)

    if inplace:
        adata.var[key_added] = roc_auc
    else:
        return roc_auc


def test_train_split(
    adata,
    *,
    test_size=None,
    train_size=None,
    random_state=0,
    shuffle: bool = True,
    stratify_col: str = None
):
    """Split an anndata object into a test and a train set.

    This mimics the sklearn API. `stratify` is a column from `adata.obs` that contains the grouping variable
    """
    if stratify_col is not None:
        stratify = adata.obs[stratify_col].values
    else:
        stratify = None
    obs_train, obs_test = sklearn.model_selection.train_test_split(
        adata.obs_names,
        test_size=test_size,
        train_size=train_size,
        shuffle=shuffle,
        random_state=random_state,
        stratify=stratify,
    )
    return adata[obs_train, :].copy(), adata[obs_test, :].copy()


def stratified_k_fold_split(
    adata, *, stratify_col: str, n_splits=5, shuffle=True, random_state=0
):
    """
    Stratified K fold split of an anndata object

    Yields
    ------
    adata_train, adata_test tuples
    """
    kf = sklearn.model_selection.StratifiedKFold(
        n_splits=n_splits, shuffle=shuffle, random_state=random_state
    )
    for train_obs, test_obs in kf.split(adata.obs_names, adata.obs[stratify_col]):
        yield adata[train_obs, :].copy(), adata[test_obs, :].copy()


def grid_search_cv(
    adata, *, replicate_col, label_col, positive_class, n_splits=5, param_grid: dict
):
    """Perform a cross-validation grid search with an MCPSignature regressor.

    Specialty is that we train on classes, but we evaluate on continuous fractions,
    both derived from a different fold of the same single cell dataset.

    Parameters
    ----------
    adata
        single-cell anndata object. Pseudobulk will be generated automatically.
        X should contain raw counts, pseudobulk will bere-normalized to log(CPM).
    replicate_col
        column in adata.obs that contains information about biological replicates (e.g. patient)
    label_col
        column in adata.obs that contains the cell-type annotation for which to generate signatures
    positive_class
        The label from label_col for which to generate signatures
    n_splits
        number of splits to use in cross-validation
    param_grid
        dictionary with lists for each parameters, e.g.

        ```python
        param_grid = [
            {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']},
        ]
        ```
    """
    replicates = np.unique(adata.obs[replicate_col])
    skf = sklearn.model_selection.KFold(n_splits=n_splits, shuffle=True, random_state=0)

    def _get_grid(param_grid):
        """Yield dictionaries {param1: v1, param2: v2, param3: v3, ...}"""
        for keys, values in zip(
            itertools.repeat(param_grid.keys()), itertools.product(*param_grid.values())
        ):
            yield {k: v for k, v in zip(keys, values)}

    grid = list(_get_grid(param_grid))

    results = []  # tuples (fold, params, score, n_genes)
    for i, (reps_train_idx, reps_test_idx) in tqdm(
        enumerate(skf.split(replicates)), total=n_splits
    ):
        reps_train, reps_test = replicates[reps_train_idx], replicates[reps_test_idx]
        print("Generating Pseudobulk")
        pb_train = pseudobulk(
            adata[adata.obs[replicate_col].isin(reps_train), :].copy(),
            groupby=[replicate_col, label_col],
        )
        pb_test = pseudobulk(
            adata[adata.obs[replicate_col].isin(reps_test), :].copy(),
            groupby=[replicate_col],  # do not include label column here
        )
        pb_test.obs.set_index("patient", inplace=True)
        pb_test.obs["true_frac"] = (
            adata.obs.groupby("patient", observed=True)
            .apply(lambda x: x[label_col].value_counts(normalize=True, dropna=False))
            .unstack()[positive_class]
        )

        # renormalize
        for ad in [pb_train, pb_test]:
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad, base=2)

        print("Computing per-gene metrics")
        pb_train = MCPSignatureRegressor.prepare_anndata(
            pb_train, obs_col=label_col, positive_class=positive_class
        )

        print("Evaluating grid")
        for params in grid:
            mcpr = MCPSignatureRegressor(**params)
            mcpr.fit(pb_train)
            try:
                y_pred = mcpr.predict(pb_test)
                score_pearson = scipy.stats.pearsonr(pb_test.obs["true_frac"], y_pred)
                n_genes = len(mcpr.signature_genes)
            except ValueError:
                score_pearson = np.nan
                n_genes = 0
            results.append(
                {
                    "fold": i,
                    **params,
                    "score_pearson": score_pearson,
                    "n_genes": n_genes,
                }
            )

    return results


class MCPSignatureRegressor:
    def __init__(self, min_fc=2, min_sfc=1.5, min_auroc=0.97):
        """
        Implements the MCP-counter appraoch for signature generation with an sklearn-like API.
        """
        self.min_fc = min_fc
        self.min_sfc = min_sfc
        self.min_auroc = min_auroc
        self._signature_genes = None

    @staticmethod
    def prepare_anndata(adata, *, obs_col: str, positive_class: str):
        """Pre-compute scores on an anndata object. Saves expensive re-computations."""
        adata = adata.copy()
        roc_auc(adata, obs_col=obs_col, positive_class=positive_class, inplace=True)
        fold_change(adata, obs_col=obs_col, positive_class=positive_class, inplace=True)
        specific_fold_change(
            adata, obs_col=obs_col, positive_class=positive_class, inplace=True
        )
        adata.uns["MCPSignatureRegressor"] = {"prepared": True}
        return adata

    def fit(self, X: AnnData):
        """
        Compute scores and select signature genes based on the specified cutoffs.

        Parameters
        ----------
        X
            anndata object. It is recommended that values in .X are log-normalized.
        obs_col
            column in adata.obs that contains the class labels
        positive_class
            label in `obs_col` that contains the positive class label
        """
        try:
            assert X.uns["MCPSignatureRegressor"]["prepared"]
        except (KeyError, AssertionError):
            raise ValueError(
                "Scores need to be pre-calculated with MCPSignatureRegressor.prepare_anndata"
            )

        self._signature_genes = X.var_names[
            (X.var["roc_auc"] >= self.min_auroc)
            & (X.var["fold_change"] >= self.min_fc)
            & (X.var["specific_fold_change"] >= self.min_sfc)
        ].tolist()

    def predict(self, X: AnnData):
        """Return the signature score for all samples in X.

        This implements an approach from the LCAM paper:
         * z-score each gene across patients
         * the mean z-score across genes is the signature score.
        """
        if self._signature_genes is None:
            raise ValueError("Model is not fitted!")
        elif not len(self._signature_genes):
            raise ValueError("Model has an empty signature!")
        X = X[:, self._signature_genes].copy()
        X.X = scipy.stats.zscore(X.X, axis=0)
        return np.mean(X.X, axis=1)

    @property
    def signature_genes(self):
        return self._signature_genes

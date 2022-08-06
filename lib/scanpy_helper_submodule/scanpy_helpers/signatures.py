"""
Functions to find and test gene signatures.

Inspired by MCPCounter and BioQC.

The original cutoffs used in the MCP counter paper are
 * fold change >= ``min_fc`` (default=2)
 * specific fold change >= ``min_sfc`` (default=1.5)
 * AUC ROC >= ``min_auc`` (default=0.97)
"""

from typing import List, Optional
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
import altair as alt
import pandas as pd


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
    """

    def _roc_auc(arr, mask):
        return roc_auc_score(mask, arr)

    positive_mask = adata.obs[obs_col] == positive_class
    roc_auc = np.apply_along_axis(_roc_auc, 0, adata.X, positive_mask)

    if inplace:
        adata.var[key_added] = roc_auc
    else:
        return roc_auc


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
    def prepare_anndata(adata, *, label_col: str, positive_class: str):
        """Pre-compute scores on an anndata object. Saves expensive re-computations.

        Parameters
        ----------
        adata
            pseudobulk anndata
        label_col
            column in adata.obs that contains the class labels
        positive_class
            label in `obs_col` that contains the positive class label
        """
        adata = adata.copy()
        roc_auc(adata, obs_col=label_col, positive_class=positive_class, inplace=True)
        fold_change(
            adata, obs_col=label_col, positive_class=positive_class, inplace=True
        )
        specific_fold_change(
            adata, obs_col=label_col, positive_class=positive_class, inplace=True
        )
        adata.uns["MCPSignatureRegressor"] = {"prepared": True}
        return adata

    def _check_adata(self, adata):
        """Check if anndata object was prepared with `prepare_anndata()`"""
        try:
            assert adata.uns[self.__class__.__name__]["prepared"]
        except (KeyError, AssertionError):
            raise ValueError(
                f"Scores need to be pre-calculated with {self.__class__.__name__}.prepare_anndata"
            )

    def fit(self, X: AnnData):
        """
        Compute scores and select signature genes based on the specified cutoffs.

        Parameters
        ----------
        X
            Pseudobulk anndata object. It is recommended that values in .X are log-normalized.
            Scores need to be precalculated with prepare_anndata().
        """
        self._check_adata(X)

        self._signature_genes = X.var.loc[
            (X.var["roc_auc"] >= self.min_auroc)
            & (X.var["fold_change"] >= self.min_fc)
            & (X.var["specific_fold_change"] >= self.min_sfc),
            :,
        ]

    def predict(self, X: AnnData):
        """Return the signature score for all samples in X.

        This implements an approach from the LCAM paper:
         * z-score each gene across patients
         * the mean z-score across genes is the signature score.
        """
        X = X[:, self.signature_genes].copy()
        X = X.X

        X = scipy.stats.zscore(X, axis=0)
        # remove genes without a z-score (e.g. all zero across all patients)
        X = X[:, ~np.any(np.isnan(X), axis=0)]

        return np.mean(X, axis=1)

    @property
    def signature_genes(self):
        """A list of the signature genes in no particular order"""
        return self._signature_genes.index.tolist()

    @property
    def signature_df(self):
        """A dataframe with the signature genes and their associated scores"""
        return self._signature_genes

    @staticmethod
    def score_pearson(true_fractions, predicted_scores):
        """Compute pearson correlaltion between true fractions and prediced cell-type scores"""
        return scipy.stats.pearsonr(true_fractions, predicted_scores)[0]

    # @staticmethod
    # def score_lm_r2(true_fractions, predicted_scores, covaraite):
    #     """Use r2 of a linear model to score the predictions. Covariate usually is a caregorical variable
    #     for confounding factors such as dataset.
    #     """
    #     pass


class HierarchicalMCPSignatureRegressor(MCPSignatureRegressor):
    def __init__(
        self,
        min_fc_coarse=2,
        min_sfc_coarse=1.5,
        min_auroc_coarse=0.97,
        min_fc_fine=1,
        min_sfc_fine=1,
        min_auroc_fine=0.8,
    ):
        """An extension to the MCP signature regressor that takes into account different cutoffs
        for two different hierarchy levels. The "coarse grained" (="low res") cutoffs are
        to distinguish major cell types. The "fine grained" (="high res") cutoffs are to distinguish
        subtypes of a cell-type.

        Marker genes will be selected for the signature that pass both the threshold for low res and
        for high res.
        """
        self._signature_genes = None
        self.min_fc_coarse = min_fc_coarse
        self.min_sfc_coarse = min_sfc_coarse
        self.min_auroc_coarse = min_auroc_coarse
        self.min_fc_fine = min_fc_fine
        self.min_sfc_fine = min_sfc_fine
        self.min_auroc_fine = min_auroc_fine

    @staticmethod
    def prepare_anndata(
        adata,
        *,
        label_col_coarse,
        label_col_fine,
        positive_class_coarse,
        positive_class_fine,
    ):
        roc_auc(
            adata,
            obs_col=label_col_coarse,
            positive_class=positive_class_coarse,
            inplace=True,
            key_added="roc_auc_coarse",
        )
        fold_change(
            adata,
            obs_col=label_col_coarse,
            positive_class=positive_class_coarse,
            inplace=True,
            key_added="fold_change_coarse",
        )
        specific_fold_change(
            adata,
            obs_col=label_col_coarse,
            positive_class=positive_class_coarse,
            inplace=True,
            key_added="specific_fold_change_coarse",
        )
        adata_fine = adata[
            adata.obs[label_col_coarse] == positive_class_coarse, :
        ].copy()
        adata.var["roc_auc_fine"] = roc_auc(
            adata_fine,
            obs_col=label_col_fine,
            positive_class=positive_class_fine,
            inplace=False,
        )
        adata.var["fold_change_fine"] = fold_change(
            adata_fine,
            obs_col=label_col_fine,
            positive_class=positive_class_fine,
            inplace=False,
        )
        adata.var["specific_fold_change_fine"] = specific_fold_change(
            adata_fine,
            obs_col=label_col_fine,
            positive_class=positive_class_fine,
            inplace=False,
        )
        adata.uns["HierarchicalMCPSignatureRegressor"] = {"prepared": True}
        return adata

    def fit(self, X: AnnData):
        self._check_adata(X)

        self._signature_genes = X.var.loc[
            (X.var["roc_auc_coarse"] >= self.min_auroc_coarse)
            & (X.var["fold_change_coarse"] >= self.min_fc_coarse)
            & (X.var["specific_fold_change_coarse"] >= self.min_sfc_coarse)
            & (X.var["roc_auc_fine"] >= self.min_auroc_fine)
            & (X.var["fold_change_fine"] >= self.min_fc_fine)
            & (X.var["specific_fold_change_fine"] >= self.min_sfc_fine),
            :,
        ]
        pass


def train_test_split(adata, *, replicate_col: str):
    """
    Split anndata object by replicates into test and train dataset.

    Parameters
    ----------
    adata
        anndata object
    replicate_col
        column in adata.obs that contains information about biological replicates (e.g. patient).
        The test/train dataset will be split by this column.
    """
    replicates = adata.obs[replicate_col].unique()
    rep_train, rep_test = sklearn.model_selection.train_test_split(
        replicates, random_state=0
    )
    return (
        adata[adata.obs[replicate_col].isin(rep_train), :].copy(),
        adata[adata.obs[replicate_col].isin(rep_test), :].copy(),
    )


def _get_grid(param_grid):
    """Yield dictionaries {param1: v1, param2: v2, param3: v3, ...}"""
    for keys, values in zip(
        itertools.repeat(param_grid.keys()), itertools.product(*param_grid.values())
    ):
        yield {k: v for k, v in zip(keys, values)}


def grid_search_cv(
    adata: AnnData,
    *,
    replicate_col: str,
    label_col: str,
    positive_class: str,
    label_col_fine: Optional[str] = None,
    positive_class_fine: Optional[str] = None,
    label_col_eval: Optional[str] = None,
    positive_class_eval: Optional[str] = None,
    n_splits: int = 5,
    model_class: type = MCPSignatureRegressor,
    param_grid: dict,
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
    label_col_fine
        Column with the fine-grained cell-type annotation labels. Only relevant for HierarchicalMCPSignatureRegressor.
    positive_class_fine
        The label from label_col_fine for which to generate signatures. Only relevant for HierarchicalMCPSignatureRegressor.
    label_col_eval
        column in adata.obs that contains the cell-type annotation which to use for the evaluation. Will
        be set to label_col if omitted.
    positive_class_eval
        The label in label_col_eval representing the positive class to be used for the evaluation. Will be set
        to positive_class if omitted.
    n_splits
        number of splits to use in cross-validation
    param_grid
        dictionary with lists for each parameters, e.g.

        ```python
        param_grid = [
            {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']},
        ]
        ```
    model_class
        The model to evaluate. Currently supported are MCPSignatureRegressor and HierarchicalMCPSignatureRegressor
    """
    if model_class not in [MCPSignatureRegressor, HierarchicalMCPSignatureRegressor]:
        raise ValueError("Unsupported model_class")

    if label_col_eval is None:
        label_col_eval = label_col
    if positive_class_eval is None:
        positive_class_eval = positive_class
    pseudobulk_groupby = (
        [label_col, label_col_fine] if label_col_fine is not None else [label_col]
    )

    replicates = np.unique(adata.obs[replicate_col])
    skf = sklearn.model_selection.KFold(n_splits=n_splits, shuffle=True, random_state=0)

    grid = list(_get_grid(param_grid))

    results = []  # tuples (fold, params, score, n_genes)
    for i, (reps_train_idx, reps_test_idx) in tqdm(
        enumerate(skf.split(replicates)), total=n_splits
    ):
        reps_train, reps_test = replicates[reps_train_idx], replicates[reps_test_idx]
        print("Generating Pseudobulk")
        pb_train = pseudobulk(
            adata[adata.obs[replicate_col].isin(reps_train), :].copy(),
            groupby=[replicate_col] + pseudobulk_groupby,
        )
        pb_test = pseudobulk(
            adata[adata.obs[replicate_col].isin(reps_test), :].copy(),
            groupby=[replicate_col],  # do not include label column here
        )
        pb_test.obs.set_index(replicate_col, inplace=True)
        pb_test.obs["true_frac"] = (
            adata.obs.groupby(replicate_col, observed=True)
            .apply(
                lambda x: x[label_col_eval].value_counts(normalize=True, dropna=False)
            )
            .unstack()[positive_class_eval]
            .fillna(0)
        )

        # renormalize
        for ad in [pb_train, pb_test]:
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad, base=2)

        print("Computing per-gene metrics")
        if model_class == MCPSignatureRegressor:
            pb_train = model_class.prepare_anndata(
                pb_train,
                label_col=label_col,
                positive_class=positive_class,
            )
        else:
            pb_train = model_class.prepare_anndata(
                pb_train,
                label_col_coarse=label_col,
                positive_class_coarse=positive_class,
                label_col_fine=label_col_fine,
                positive_class_fine=positive_class_fine,
            )

        print("Evaluating grid")
        for params in grid:
            mcpr = model_class(**params)
            mcpr.fit(pb_train)
            n_genes = len(mcpr.signature_genes)
            if not n_genes:
                score_pearson = np.nan
            else:
                y_pred = mcpr.predict(pb_test)
                score_pearson = mcpr.score_pearson(
                    pb_test.obs["true_frac"].values, y_pred
                )

            results.append(
                {
                    "fold": i,
                    **params,
                    "score_pearson": score_pearson,
                    "n_genes": n_genes,
                }
            )

    return results


def refit_and_evaluate(
    model: MCPSignatureRegressor,
    adata_train: AnnData,
    adata_test: AnnData,
    *,
    replicate_col: str,
    label_col: str,
    covariate_col: str,
    positive_class: str,
):
    """
    Refit a model after cross validation using the optimal parameters and evaluate it on the independent test
    dataset.

    Parameters
    ----------
    model
        an instance of a MCPSignatureRegressor model
    adata_train
        single-cell anndata object with the training data
    adata_test
        single-cell anndata object with the test data
    replicate_col
        column in adata.obs that contains information about biological replicates (e.g. patient)
    label_col
        column in adata.obs that contains the cell-type annotation for which to generate signatures
    positive_class
        label in `label_col` that contains the positive class label

    Returns
    -------
    scores
        The final scores computed on the independent test dataset
    """
    pb_train = pseudobulk(adata_train, groupby=[replicate_col, label_col])
    pb_test = pseudobulk(adata_test, groupby=[replicate_col, covariate_col])
    pb_test.obs.set_index(replicate_col, inplace=True)
    pb_test.obs["true_frac"] = (
        adata_test.obs.groupby(replicate_col, observed=True)
        .apply(lambda x: x[label_col].value_counts(normalize=True, dropna=False))
        .unstack()[positive_class]
        .fillna(0)
    )

    # renormalize
    for ad in [pb_train, pb_test]:
        sc.pp.normalize_total(ad, target_sum=1e6)
        sc.pp.log1p(ad, base=2)

    pb_train = model.prepare_anndata(
        pb_train, label_col=label_col, positive_class=positive_class
    )

    model.fit(pb_train)

    y_pred = model.predict(pb_test)

    result_df = pd.DataFrame().assign(
        true_frac=pb_test.obs["true_frac"],
        predicted_score=y_pred,
        dataset=pb_test.obs["dataset"].values,
    )

    score_pearson = model.score_pearson(pb_test.obs["true_frac"].values, y_pred)

    return {
        "n_genes": len(model.signature_genes),
        "score_pearson": score_pearson,
        "result_df": result_df,
    }

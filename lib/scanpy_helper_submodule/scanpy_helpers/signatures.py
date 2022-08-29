"""
Functions to find and test gene signatures.

Inspired by MCPCounter and BioQC.

The original cutoffs used in the MCP counter paper are
 * fold change >= ``min_fc`` (default=2)
 * specific fold change >= ``min_sfc`` (default=1.5)
 * AUC ROC >= ``min_auc`` (default=0.97)
"""

from typing import Callable, Dict, List, Mapping, Optional, Sequence
import numpy as np
from sklearn.metrics import roc_auc_score
import sklearn.model_selection
from anndata import AnnData
import itertools
import scipy.stats
from .pseudobulk import pseudobulk
import scanpy as sc
import altair as alt
import pandas as pd
from tqdm.contrib.concurrent import process_map


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
    def __init__(
        self,
        min_fc: float = 2,
        min_sfc: float = 1.5,
        min_auroc: float = 0.97,
        constraint: Optional[Callable[[pd.DataFrame], pd.DataFrame]] = None,
    ):
        """
        Implements the MCP-counter appraoch for signature generation with an sklearn-like API.

        min_fc
            minimum fold change of a gene to count as marker gene
        min_sfc
            minimum specific fold change of a gene to count as marker gene
        min_auroc
            minimum auroc of a gene to count as marker gene
        constraint
            An additional constraint on signature genes. This is a function to be executed on `adata.var`.
            Takes a dataframe as input and returns a filtered data frame with the same columns.
        """
        self.min_fc = min_fc
        self.min_sfc = min_sfc
        self.min_auroc = min_auroc
        self._signature_genes = None
        self.constraint = constraint

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

        tmp_signature_genes = X.var.loc[
            (X.var["roc_auc"] >= self.min_auroc)
            & (X.var["fold_change"] >= self.min_fc)
            & (X.var["specific_fold_change"] >= self.min_sfc),
            :,
        ]

        if self.constraint is not None:
            tmp_signature_genes = self.constraint(tmp_signature_genes)

        self._signature_genes = tmp_signature_genes

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


def _grid_search_cv_execute_fold(
    i,
    reps_train_labels,
    reps_test_labels,
    pb_train,
    pb_test,
    replicate_col,
    label_col,
    positive_class,
    grid,
):
    """Execute a single fold for the grid search. Used for parallelization"""
    results = []  # tuples (fold, params, score, n_genes)

    pb_train = pb_train[pb_train.obs[replicate_col].isin(reps_train_labels), :]
    pb_test = pb_test[pb_test.obs_names.isin(reps_test_labels), :]

    print("Computing per-gene metrics")
    pb_train = MCPSignatureRegressor.prepare_anndata(
        pb_train,
        label_col=label_col,
        positive_class=positive_class,
    )

    print("Evaluating grid")
    for params in grid:
        mcpr = MCPSignatureRegressor(**params)
        mcpr.fit(pb_train)
        n_genes = len(mcpr.signature_genes)
        if not n_genes:
            score_pearson = np.nan
        else:
            y_pred = mcpr.predict(pb_test)
            score_pearson = mcpr.score_pearson(pb_test.obs["true_frac"].values, y_pred)

        results.append(
            {
                "fold": i,
                **params,
                "score_pearson": score_pearson,
                "n_genes": n_genes,
            }
        )

    return results


def grid_search_cv(
    adata: AnnData,
    *,
    replicate_col: str,
    label_col: str,
    positive_class: str,
    n_splits: int = 5,
    param_grid: dict,
    raw_results: bool = False,
    n_jobs=None,
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
    raw_results
        If True, return a dictionary with the results of each fold. Otherwise aggregate the results into a dataframe
        with one row per parameter. The dataframe will contain this function's parameters in `attrs`
        which is useful for posteriority and for downstream functions.
    """
    replicates = np.unique(adata.obs[replicate_col])
    skf = sklearn.model_selection.KFold(n_splits=n_splits, shuffle=True, random_state=0)
    reps_train_labels, reps_test_labels = zip(
        *[
            (replicates[train_idx], replicates[test_idx])
            for train_idx, test_idx in skf.split(replicates)
        ]
    )

    grid = list(_get_grid(param_grid))

    # pb train, test are both grouped by patient (they will be split in test/train patients later)
    # The train pseudobulk is additionally split into the cell types, while the test pseudobulk
    # has all cell-types mixed bet has the original cell-type fractions annotated for validation.
    print("Generating Pseudobulk")
    pb_train = pseudobulk(
        adata,
        groupby=[replicate_col, label_col],
    )
    pb_test = pseudobulk(
        adata,
        groupby=[replicate_col],  # do not include label column here
    )
    pb_test.obs.set_index(replicate_col, inplace=True)
    pb_test.obs["true_frac"] = (
        adata.obs.groupby(replicate_col, observed=True)
        .apply(lambda x: x[label_col].value_counts(normalize=True, dropna=False))
        .unstack()[positive_class]
        .fillna(0)
    )

    # renormalize
    for ad in [pb_train, pb_test]:
        sc.pp.normalize_total(ad, target_sum=1e6)
        sc.pp.log1p(ad, base=2)

    all_results = process_map(
        _grid_search_cv_execute_fold,
        list(enumerate(reps_train_labels)),
        reps_train_labels,
        reps_test_labels,
        itertools.repeat(pb_train),
        itertools.repeat(pb_test),
        itertools.repeat(replicate_col),
        itertools.repeat(label_col),
        itertools.repeat(positive_class),
        itertools.repeat(grid),
        max_workers=n_jobs,
    )
    all_results = list(itertools.chain.from_iterable(all_results))

    if raw_results:
        return all_results
    else:
        return results_to_df(
            all_results,
            attrs={
                "replicate_col": replicate_col,
                "label_col": label_col,
                "positive_class": positive_class,
            },
        )


def results_to_df(results: List[Dict], attrs: Optional[Mapping] = None):
    """Make result dataframe from dictionaries. Optionally stores attributes in df.attrs."""
    df = pd.DataFrame(results).drop(columns="fold")

    groupby_cols = [x for x in df.columns if x not in ["n_genes", "score_pearson"]]

    df = (
        df.groupby(groupby_cols)
        .agg(lambda x: x.mean(skipna=False))
        .sort_values("score_pearson", ascending=False)
        .reset_index()
    )
    if attrs is not None:
        for k, v in attrs.items():
            df.attrs[k] = v

    return df


def refit_and_evaluate(
    model: MCPSignatureRegressor,
    adata_train: AnnData,
    adata_test: AnnData,
    *,
    replicate_col: str,
    label_col: str,
    covariate_cols: Sequence[str] = (),
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
    covariate_cols
        Additional covariates to include into pseudobulk groupby

    Returns
    -------
    scores
        The final scores computed on the independent test dataset
    """
    pb_train = pseudobulk(adata_train, groupby=[replicate_col, label_col])
    pb_test = pseudobulk(adata_test, groupby=[replicate_col] + covariate_cols)
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


def refit_evaluate_plot(
    scores, *, i=0, adata_train, adata_test, plot_covariate="dataset"
):
    """A high-level wrapper around refit_and_evaluate"""
    params = {
        k: v for k, v in scores.iloc[i].items() if k not in ["n_genes", "score_pearson"]
    }
    print(params)
    mcpr = MCPSignatureRegressor(**params)
    res = refit_and_evaluate(
        mcpr,
        adata_train,
        adata_test,
        replicate_col=scores.attrs["replicate_col"],
        label_col=scores.attrs["label_col"],
        positive_class=scores.attrs["positive_class"],
        covariate_cols=[plot_covariate],
    )
    print(f"n_genes: {res['n_genes']}, score_pearson: {res['score_pearson']}")
    alt.Chart(res["result_df"]).mark_circle().encode(
        x="true_frac", y="predicted_score", color=plot_covariate
    ).display()
    return mcpr


def plot_metric_strip(
    pb, markers, *, metric_key="auroc", metric_title="AUROC", domain=[0.5, 1], top=10
):
    """A one-row heatmap accompaniying the heatmap generated by plot_markers showing the marker
    quality."""
    metric = []
    genes = []
    cts = []
    for ct, tmp_markers in markers.items():
        for gene in tmp_markers[:top]:
            metric.append(pb.var.loc[gene, f"{ct}_{metric_key}"])
            genes.append(gene)
            cts.append(ct)

    tmp_df = pd.DataFrame().assign(marker=genes, cell_type=cts, metric=metric)

    return (
        alt.Chart(tmp_df)
        .mark_rect()
        .encode(
            x=alt.X("marker", sort=None),
            color=alt.Color(
                "metric",
                scale=alt.Scale(scheme="viridis", domain=domain),
                title=metric_title,
            ),
        )
    )


def plot_markers(pb, groupby, markers, *, top=10, **kwargs):
    return sc.pl.matrixplot(
        pb,
        groupby=groupby,
        var_names={k: v[:top] for k, v in markers.items()},
        cmap="bwr",
        **kwargs,
    )

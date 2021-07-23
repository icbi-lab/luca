"""Visualize scanpy QC metrics"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd


def get_stats_df(adata, dataset_id):
    """
    Generate a dataframe with QC stats
    """
    return pd.DataFrame().assign(
        dataset_id=[dataset_id],
        min_genes=[np.min(adata.obs["n_genes_by_counts"])],
        max_genes=[np.max(adata.obs["n_genes_by_counts"])],
        min_counts=[np.min(adata.obs["total_counts"])],
        max_counts=[np.max(adata.obs["total_counts"])],
        min_pct_mito=[np.min(adata.obs["pct_counts_mito"])],
        max_pct_mito=[np.max(adata.obs["pct_counts_mito"])],
        n_obs=len(adata.obs_names),
        n_var=len(adata.var_names),
    )


def plot_qc_metrics(
    adata,
    cumulative=False,
    min_genes=None,
    max_genes=None,
    min_counts=None,
    max_pct_mito=None,
    max_counts=None,
    show=True,
):
    """
    Plots the QC metrics generated with scanpy.pp.calculate_qc_metrics.

    calculate_qc_metrics needs to be ran with `qc_vars=("mito", )`.

    The min/max cutoffs are only used for plotting cutoff lines.
    """
    KEYS = (
        "total_counts",
        "total_counts (< 5000)",
        "n_genes_by_counts (< 1000)",
        "n_genes_by_counts",
        "pct_counts_mito",
    )
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(17, 10))
    axs = {k: ax for k, ax in zip(KEYS, axs.flatten())}

    def _add_line(ax, val):
        if val is not None:
            if cumulative:
                ax.hlines(
                    y=val, color="red", linewidth=1, xmin=0, xmax=ax.get_xlim()[1]
                )
            else:
                ax.vlines(
                    x=val, color="red", linewidth=1, ymin=0, ymax=ax.get_ylim()[1]
                )

    for key, ax in axs.items():
        if key == "total_counts (< 5000)":
            values = adata.obs["total_counts"][adata.obs["total_counts"] < 5000]
        elif key == "n_genes_by_counts (< 1000)":
            values = adata.obs["n_genes_by_counts"][
                adata.obs["n_genes_by_counts"] < 1000
            ]
        else:
            values = adata.obs[key]

        if not cumulative:
            sns.histplot(values, bins=100, ax=ax)
        else:
            ranks = values.rank(ascending=True, method="first")
            ax.scatter(x=ranks, y=values, marker=".", s=2)

        ax.set_xlabel(key)

    if cumulative:
        axs["total_counts"].set_yscale("log")
        axs["total_counts (< 5000)"].set_yscale("log")

    _add_line(axs["total_counts"], min_counts)
    _add_line(axs["total_counts (< 5000)"], min_counts)
    _add_line(axs["total_counts"], max_counts)
    _add_line(axs["n_genes_by_counts"], min_genes)
    _add_line(axs["n_genes_by_counts (< 1000)"], min_genes)
    _add_line(axs["n_genes_by_counts"], max_genes)
    _add_line(axs["pct_counts_mito"], max_pct_mito)

    if show:
        plt.show()
    else:
        return fig

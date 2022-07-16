"""Plotting functions for group comparisons"""

import altair as alt
import pandas as pd
import numpy as np


def plot_lm_result_altair(
    df,
    p_cutoff=0.1,
    p_col="fdr",
    x="variable",
    y="group",
    color="coef",
    title="heatmap",
    cluster=False,
    value_max=None,
    configure=lambda x: x.configure_mark(opacity=1),
    cmap="redblue",
    reverse=True,
    domain=lambda x: [-x, x],
    order=None,
):
    """
    Plot a results data frame of a comparison as a heatmap
    """
    df_filtered = df.loc[lambda _: _[p_col] < p_cutoff, :]
    df_subset = df.loc[
        lambda _: _[x].isin(df_filtered[x].unique()) & _[y].isin(df[y].unique())
    ]
    if not df_subset.shape[0]:
        print("No values to plot")
        return

    if order is None:
        order = "ascending"
        if cluster:
            from scipy.cluster.hierarchy import linkage, leaves_list

            values_df = df_subset.pivot(index=y, columns=x, values=color)
            order = values_df.columns.values[
                leaves_list(
                    linkage(values_df.values.T, method="average", metric="euclidean")
                )
            ]

    def _get_significance(fdr):
        if fdr < 0.001:
            return "< 0.001"
        elif fdr < 0.01:
            return "< 0.01"
        elif fdr < 0.1:
            return "< 0.1"
        else:
            return np.nan

    df_subset["FDR"] = pd.Categorical([_get_significance(x) for x in df_subset[p_col]])

    if value_max is None:
        value_max = max(
            abs(np.nanmin(df_subset[color])), abs(np.nanmax(df_subset[color]))
        )
    # just setting the domain in altair will lead to "black" fields. Therefore, we constrain the values themselves.
    df_subset[color] = np.clip(df_subset[color], *domain(value_max))
    return configure(
        alt.Chart(df_subset, title=title)
        .mark_rect()
        .encode(
            x=alt.X(x, sort=order),
            y=y,
            color=alt.Color(
                color,
                scale=alt.Scale(scheme=cmap, reverse=reverse, domain=domain(value_max)),
            ),
        )
        + alt.Chart(df_subset.loc[lambda x: ~x["FDR"].isnull()])
        .mark_point(color="white", filled=True, stroke="black", strokeWidth=0)
        .encode(
            x=alt.X(x, sort=order),
            y=y,
            size=alt.Size(
                "FDR:N",
                scale=alt.Scale(
                    domain=["< 0.001", "< 0.01", "< 0.1"],
                    range=4 ** np.array([3, 2, 1]),
                ),
            ),
        )
    )

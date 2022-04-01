from typing import List, Optional
import pandas as pd
import numpy as np
import scanpy as sc
import importlib.resources as pkg_resources
from . import assets
from tqdm import tqdm
from scanpy import logging
import matplotlib.pyplot as plt


class AnnotationHelper:
    def __init__(self, sheet="lung", markers=None):
        """sheet is not used anymore, uses the internal CSV file as default now."""
        if markers is not None:
            self.markers = markers
        else:
            self.markers = pd.read_csv(
                pkg_resources.open_text(assets, "lung_markers.csv")
            )
            # self.markers = pd.read_csv(
            #     f"https://docs.google.com/spreadsheets/d/1beW-9oeM31P50NFvNLVvsdXlfh_tOjsmIrnwV2ZlxDU/gviz/tq?tqx=out:csv&sheet={sheet}"
            # )

    def get_markers(self, filter_cell_type: Optional[List] = None):
        if filter_cell_type is None:
            return self.markers.copy()
        else:
            mask = np.zeros(self.markers.shape[0])
            for filter in filter_cell_type:
                mask = mask | self.markers["cell_type"].str.lower().str.contains(
                    filter.lower()
                )
            return self.markers.loc[mask, :]

    @staticmethod
    def get_marker_dict(markers, adata):
        marker_dict = dict()
        for cell_type, gene in zip(markers["cell_type"], markers["gene_identifier"]):
            var_names = adata.var_names if adata.raw is None else adata.raw.var_names
            if gene not in var_names:
                print(f"Marker {gene} for {cell_type} not in adata.var")
                continue
            try:
                marker_dict[cell_type].append(gene)
            except KeyError:
                marker_dict[cell_type] = [gene]

        return marker_dict

    def plot_dotplot(
        self,
        adata,
        *,
        markers=None,
        filter_cell_type: Optional[List] = None,
        groupby="leiden",
    ):
        if markers is not None:
            tmp_markers = markers
        else:
            tmp_markers = self.get_markers(filter_cell_type)
        marker_dict = self.get_marker_dict(tmp_markers, adata)

        sc.pl.dotplot(adata, var_names=marker_dict, groupby=groupby)

    def plot_umap(
        self, adata, *, markers=None, filter_cell_type: Optional[List] = None, **kwargs
    ):
        """Make a UMAP plot for each marker. Filter markers by explicitly
        specifying them or providing a filter keyword for the marker table.
        Kwargs are passed to sc.pl.umap()"""
        if markers is not None:
            tmp_markers = markers
        else:
            tmp_markers = self.get_markers(filter_cell_type)

        for ct, group in tmp_markers.groupby("cell_type"):
            print(ct)
            var_names = adata.var_names if adata.raw is None else adata.raw.var_names
            genes = set(var_names) & set(group["gene_identifier"].values)
            sc.pl.umap(adata, color=genes, **kwargs)

    def score_cell_types(
        self,
        adata,
        *,
        markers=None,
        filter_cell_type: Optional[List] = None,
        prefix="ct_",
        **kwargs,
    ):
        tmp_markers = self.get_markers(filter_cell_type)
        for ct, group in tqdm(tmp_markers.groupby("cell_type")):
            var_names = adata.var_names if adata.raw is None else adata.raw.var_names
            genes = list(set(var_names) & set(group["gene_identifier"].values))
            sc.tl.score_genes(adata, genes, score_name=prefix + ct, **kwargs)

    def plot_dotplot_scores(
        self,
        adata,
        *,
        markers=None,
        filter_cell_type: Optional[List] = None,
        prefix="ct_",
        groupby="leiden",
        **kwargs,
    ):
        if markers is not None:
            tmp_markers = markers
        else:
            tmp_markers = self.get_markers(filter_cell_type)

        var_names = [prefix + x for x in tmp_markers["cell_type"].unique()]
        if set(var_names) - set(adata.obs.columns) != set():
            logging.info(
                "Scores not found in adata.obs. Computing scores with default parameters. "
            )
            self.score_cell_types(
                adata,
                markers=tmp_markers,
                filter_cell_type=filter_cell_type,
                prefix=prefix,
            )

        sc.pl.dotplot(
            adata,
            groupby=groupby,
            var_names=var_names,
            **kwargs,
        )

    def plot_umap_scores(
        self,
        adata,
        *,
        markers=None,
        filter_cell_type: Optional[List] = None,
        prefix="ct_",
        **kwargs,
    ):
        if markers is not None:
            tmp_markers = markers
        else:
            tmp_markers = self.get_markers(filter_cell_type)

        var_names = [prefix + x for x in tmp_markers["cell_type"].unique()]
        if set(var_names) - set(adata.obs.columns) != set():
            logging.info(
                "Scores not found in adata.obs. Computing scores with default parameters. "
            )
            self.score_cell_types(
                adata,
                markers=tmp_markers,
                filter_cell_type=filter_cell_type,
                prefix=prefix,
            )

        sc.pl.umap(
            adata,
            color=var_names,
            **kwargs,
        )

    @staticmethod
    def reprocess_adata_subset(
        adata, *, n_top_genes=2000, batch_correct=None, batch_key="batch"
    ):
        """Recompute HVG, PCA and UMAP on a cell-type cluster to perform subclustering"""
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        if batch_correct == "combat":
            sc.pp.combat(adata, key=batch_key)
        sc.tl.pca(adata)
        if batch_correct == "harmony":
            sc.external.pp.harmony_integrate(adata, key=batch_key)
            sc.pp.neighbors(adata, use_rep="X_pca_harmony")
        elif batch_correct == "bbknn":
            sc.external.pp.bbknn(adata, batch_key=batch_key)
        else:
            sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata)

    @staticmethod
    def reprocess_adata_subset_scvi(
        adata,
        *,
        n_neighbors=10,
        leiden_res=1,
        use_rep="X_scVI",
    ):
        """Recompute UMAP and leiden on a adata subset when scVI is used (no additional
        batch correction)"""
        sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=leiden_res)

    @staticmethod
    def reprocess_adata_subset_cnv(adata):
        """Recomput UMAP and leiden on a data subset with infercnvpy."""
        import infercnvpy as cnv

        cnv.tl.pca(adata)
        cnv.pp.neighbors(adata)
        cnv.tl.leiden(adata)
        cnv.tl.umap(adata)

    @staticmethod
    def annotate_cell_types(
        adata,
        cell_type_map,
        *,
        key_added="cell_type",
        default_cell_type="other",
        column="leiden",
    ):
        """Generate a cell-type column from a Mapping cell_type -> [list of clusters]"""
        res = np.full(adata.shape[0], default_cell_type, dtype=object)
        for ct, clusters in cell_type_map.items():
            clusters = [str(x) for x in clusters]
            res[adata.obs[column].isin(clusters)] = ct

        adata.obs[key_added] = res
        sc.pl.umap(adata, color=key_added)

    @staticmethod
    def integrate_back(
        adata,
        adata_subset,
        variable="cell_type",
        plt_context={"figure.figsize": (5, 5)},
    ):
        """Merge cell type annotations performed on a subset back into the main
        AnnData object"""
        adata.obs[variable] = adata.obs[variable].astype("str")
        adata.obs.loc[adata_subset.obs.index, variable] = adata_subset.obs[
            variable
        ].astype("str")
        with plt.rc_context(plt_context):
            sc.pl.umap(adata, color=variable)


def classify_cell_types_nearest_neighbors(
    adata, obs_key, *, mask_reference, mask_query, key_added, transitive=1
):
    """
    Simple cell-type classifier based on the cell 2 cell neighborhood graph.
    Performs a simple weighted majority voting based on the connectivities.

    I implemented this because the scANVI predictions performed rather poorly in my case.

    Parameters
    ----------
    adata
        merged anndata with both query and reference cells and
        a joint neighborhood graph (e.g. generated with scANVI)
    obs_key
        kye in adata.obs that contains the reference cell-types annotations
    mask_reference
        boolean array indicating reference cells
    maks_query
        boolean array indicating query cells
    key_added
        The predicted cell-type annotation is stored here. This may be the
        same as `obs_key` to add the predictions into the original column.
    transitive
        If > 0 also neighbors of neighbors are considered. The number specifies
        the steps in the graph. Using neighbors of neighbors increases the robustness.
    """
    conn = adata.obsp["connectivities"]
    for _ in range(transitive):
        conn = conn.dot(conn)
    conn_subset = conn[mask_reference, :][:, mask_query]
    ct_dummies = pd.get_dummies(adata.obs.loc[mask_reference][obs_key])
    ct_weights = conn_subset.T.dot(ct_dummies.values)
    predictions = np.array(ct_dummies.columns[np.argmax(ct_weights, axis=1)])
    adata.obs.loc[mask_query, key_added] = predictions

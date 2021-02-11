import pandas as pd
import numpy as np
import scanpy as sc


class AnnotationHelper:
    def __init__(self):
        self.markers = pd.read_csv(
            "https://docs.google.com/spreadsheets/d/1beW-9oeM31P50NFvNLVvsdXlfh_tOjsmIrnwV2ZlxDU/gviz/tq?tqx=out:csv&sheet=lung"
        )

    def get_markers(self, filter_cell_type: list = None):
        if filter_cell_type is None:
            return self.markers.copy()
        else:
            mask = np.zeros(self.markers.shape[0])
            for filter in filter_cell_type:
                mask = mask | self.markers["cell_type"].str.contains(filter)
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

    def plot_dotplot(self, adata, *, markers=None, filter_cell_type: list = None):
        if markers is not None:
            tmp_markers = markers
        else:
            tmp_markers = self.get_markers(filter_cell_type)
        marker_dict = self.get_marker_dict(tmp_markers, adata)

        sc.pl.dotplot(adata, var_names=marker_dict, groupby="leiden")

    def plot_umap(self, adata, *, markers=None, filter_cell_type: list = None):
        if markers is not None:
            tmp_markers = markers
        else:
            tmp_markers = self.get_markers(filter_cell_type)

        for ct, group in tmp_markers.groupby("cell_type"):
            print(ct)
            var_names = adata.var_names if adata.raw is None else adata.raw.var_names
            genes = set(var_names) & set(group["gene_identifier"].values)
            sc.pl.umap(adata, color=genes)

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
        """Recomput UMAP and leiden on a data subset with infercnvpy. """
        import infercnvpy as cnv

        cnv.tl.pca(adata)
        cnv.pp.neighbors(adata)
        cnv.tl.leiden(adata)
        cnv.tl.umap(adata)

    @staticmethod
    def annotate_cell_types(
        adata, cell_type_map, key_added="cell_type", default_cell_type="other"
    ):
        """Generate a cell-type column from a Mapping cell_type -> [list of clusters]"""
        adata.obs[key_added] = default_cell_type
        for ct, clusters in cell_type_map.items():
            clusters = [str(x) for x in clusters]
            adata.obs.loc[adata.obs["leiden"].isin(clusters), "cell_type"] = ct
        sc.pl.umap(adata, color="cell_type")

    @staticmethod
    def integrate_back(adata, adata_subset):
        """Merge cell type annotations performed on a subset back into the main
        AnnData object"""
        adata.obs["cell_type"] = adata.obs["cell_type"].astype("str")
        adata.obs.loc[adata_subset.obs.index, "cell_type"] = adata_subset.obs[
            "cell_type"
        ].astype("str")
        sc.pl.umap(adata, color="cell_type")

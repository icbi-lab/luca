
def modules = params.annotate_datasets.clone()

include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_COARSE }  from "../modules/local/jupyternotebook/main.nf" addParams (
    options: modules["ANNOTATE_CELL_TYPES_COARSE"]
)
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_CELL_TYPES } from "../subworkflows/neighbors_leiden_umap/main.nf" addParams(
    options: modules["NEIGHBORS_LEIDEN_UMAP_CELL_TYPES"]
)
include { PREPARE_ANNDATA as PREPARE_ANNDATA_DE_EPI }  from "../modules/local/scde/main.nf" addParams(
    options: modules["PREPARE_ANNDATA_DE_EPI"]
)
include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf" addParams(
    options: modules["SPLIT_ANNDATA"]
)
include { MAKE_PSEUDOBULK as MAKE_PSEUDOBULK_EPI }  from "../modules/local/scde/main.nf" addParams(
    options: modules["MAKE_PSEUDOBULK_EPI"]
)
include { DE_EDGER as DE_EDGER_EPI } from "../modules/local/scde/main.nf" addParams(
    options: modules["DE_EDGER_EPI"]
)
include { DE_EDGER as DE_EDGER_EPI_N_CELLS } from "../modules/local/scde/main.nf" addParams(
    options: modules["DE_EDGER_EPI_N_CELLS"]
)

/**
 * Annotate cell-types of the lung cancer atlas.
 *   - perform a coarse-grained annotation
 *   - perform DE analysis on some clusters
 *   - manually annotate sub-clusters to obtain a fine-grained cell-type annotation.
 */
workflow annotate_dataset {
    take:
        adata_integrated

    main:
    ANNOTATE_CELL_TYPES_COARSE(
        Channel.value([
            [id: "27_annotate_cell_types"],
            file("${baseDir}/analyses/20_integrate_scrnaseq_data/27_annotate_cell_types_coarse.py")
        ]),
        [:],
        adata_integrated
    )
    SPLIT_ANNDATA(
        ANNOTATE_CELL_TYPES_COARSE.out.artifacts.map{ it -> [it.baseName, it]},
        "cell_type"
    )

    NEIGHBORS_LEIDEN_UMAP_CELL_TYPES(
        SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it] },
        "X_scANVI",
        Channel.from(0.5, 0.75, 1.0, 1.5)
    )

    ch_epithelial = NEIGHBORS_LEIDEN_UMAP_CELL_TYPES.out.adata.filter{
        id, adata -> id.contains("epithelial")
    }
    PREPARE_ANNDATA_DE_EPI(
        ch_epithelial,
        "X",
        "leiden_0.50",
        ["all", "rest"]
    )
    MAKE_PSEUDOBULK_EPI(
        PREPARE_ANNDATA_DE_EPI.out.adata.flatMap{ id, adatas -> adatas }.map{ it -> [it.baseName, it]},
        "patient",
        "leiden_0.50",
        10
    )
    DE_EDGER_EPI(
        MAKE_PSEUDOBULK_EPI.out.pseudobulk,
        "leiden_0.50",
        ""
    )
    DE_EDGER_EPI_N_CELLS(
        MAKE_PSEUDOBULK_EPI.out.pseudobulk,
        "leiden_0.50",
        " + n_cells"
    )


    // emit:
    //     adata_annotated_by_cell_type = ch_adata_annotated_by_cell_type
    //     adata_annotated_cell_type_coarse = ch_adata_annotated
}

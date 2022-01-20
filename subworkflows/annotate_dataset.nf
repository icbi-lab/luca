
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_COARSE }  from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_FINE }  from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_EPI }  from "../modules/local/jupyternotebook/main.nf"
include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_CELL_TYPES } from "./neighbors_leiden_umap.nf"
include { JUPYTERNOTEBOOK as EXPORT_ATLAS }  from "../modules/local/jupyternotebook/main.nf"

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
            [id: "annotate_cell_types_coarse"],
            file("${baseDir}/analyses/30_annotate_scrnaseq_data/31_annotate_cell_types_coarse.py")
        ]),
        [:],
        adata_integrated
    )
    ch_adata_annotated = ANNOTATE_CELL_TYPES_COARSE.out.artifacts
    SPLIT_ANNDATA(
        ch_adata_annotated.map{ it -> [it.baseName, it]},
        "cell_type"
    )
    NEIGHBORS_LEIDEN_UMAP_CELL_TYPES(
        SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it] },
        "X_scANVI",
        Channel.from(0.5, 0.75, 1.0, 1.5)
    )
    ANNOTATE_CELL_TYPES_FINE(
        Channel.value([
            [id: "annotate_cell_types_fine"],
            file("${baseDir}/analyses/30_annotate_scrnaseq_data/32_annotate_cell_types_fine.py")
        ]),
        [
            "input_dir": '.',
            "main_adata": 'adata_cell_type_coarse.h5ad'
        ],
        NEIGHBORS_LEIDEN_UMAP_CELL_TYPES.out.adata.map{ id, adata -> adata }.mix(
            ch_adata_annotated
        ).collect()
    )
    ANNOTATE_CELL_TYPES_EPI(
        Channel.value([
            [id: "annotate_cell_types_epi"],
            file("${baseDir}/analyses/30_annotate_scrnaseq_data/33_epithelial_cells.py")
        ]),
        [
            "input_adata": 'adata_cell_type_coarse_epithelial_cell.umap_leiden.h5ad',
        ],
        NEIGHBORS_LEIDEN_UMAP_CELL_TYPES.out.adata.map{ id, adata -> adata }.filter(
            it -> it.name.equals("adata_cell_type_coarse_epithelial_cell.umap_leiden.h5ad")
        )
    )

    EXPORT_ATLAS(
        Channel.value([
            [id: "export_atlas"],
            file("${baseDir}/analyses/30_annotate_scrnaseq_data/35_export_atlas.py")
        ]),
        [
            "adata_annotated_fine": "adata_annotated_fine.h5ad",
            "adata_epi": "adata_epithelial.h5ad",
            "platform_metadata": "sequencing_platforms.csv",
            "patient_metadata": "patient_metadata_corrected.xlsx"
        ],
        ANNOTATE_CELL_TYPES_FINE.out.artifacts.mix(
            ANNOTATE_CELL_TYPES_EPI.out.artifacts
        ).mix(
            Channel.fromPath("$baseDir/tables/additional_patient_metadata/patient_metadata_corrected.xlsx")
        ).mix(
            Channel.fromPath("$baseDir/tables/additional_patient_metadata/sequencing_platforms.csv")
        ).collect()
    )
    ch_atlas = EXPORT_ATLAS.out.artifacts.flatten().filter{ it -> it.baseName.equals("full_atlas_annotated") }

    emit:
        final_atlas = ch_atlas
}

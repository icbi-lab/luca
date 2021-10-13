#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2
assert params.input: "Input samplesheet not specified!"

include { integrate_datasets } from "./workflows/integrate_datasets.nf"
include { annotate_dataset } from "./workflows/annotate_dataset.nf"

// include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_FINE }  from "./modules/local/jupyternotebook/main.nf" addParams (
//     options: modules["ANNOTATE_CELL_TYPES_FINE"]
// )
// include { JUPYTERNOTEBOOK as PREPARE_CELLXGENE }  from "./modules/local/jupyternotebook/main.nf" addParams (
//     options: modules["PREPARE_CELLXGENE"]
// )
// include { H5AD_TO_SEURAT }  from "./modules/local/scconversion/main.nf" addParams(
//     options: modules["H5AD_TO_SEURAT"]
// )
// include { SPLIT_ANNDATA }  from "./modules/local/scconversion/main.nf" addParams(
//     options: modules["SPLIT_ANNDATA"]
// )


workflow {

    integrate_datasets()
    annotate_dataset(integrate_datasets.out.adata_integrated)


    /// PREPARE FINAL OUTPUT

    // ANNOTATE_CELL_TYPES_FINE(
    //     Channel.value([
    //         [id: "29_annotate_cell_types_fine"],
    //         file("${baseDir}/analyses/20_integrate_scrnaseq_data/29_annotate_cell_types_fine.py")
    //     ]),
    //     [
    //         "input_dir": '.',
    //         "main_adata": 'adata_cell_type_coarse.h5ad'
    //     ],
    //     NEIGHBORS_LEIDEN_UMAP_CELL_TYPES.out.adata.map{ id, adata -> adata }.mix(
    //         ch_adata_annotated
    //     ).collect()
    // )

    // PREPARE_CELLXGENE(
    //     Channel.value([
    //         [id: "zz_prepare_cellxgene"],
    //         file("${baseDir}/analyses/zz_cellxgene/stats_and_cellxgene.py")
    //     ]),
    //     ["adata_in": "adata_annotated_fine.h5ad"],
    //     ANNOTATE_CELL_TYPES_FINE.out.artifacts
    // )

    // ch_adata_annotated_fine = ANNOTATE_CELL_TYPES_FINE.out.artifacts.flatten().filter(
    //     it -> it.name.contains("h5ad")
    // ).map{ it -> [it.baseName, it]}

    // SPLIT_ANNDATA(ch_adata_annotated_fine, "dataset")
    // H5AD_TO_SEURAT(ch_adata_annotated_fine)


}


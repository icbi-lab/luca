include { SCQC } from "../modules/local/scqc/main"
include { SCQC_MERGE_STATS } from "../modules/local/scqc_merge_stats/main.nf"
include { JUPYTERNOTEBOOK as INTEGRATE_INTO_ATLAS } from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as UPDATE_ANNOTATION } from "../modules/local/jupyternotebook/main.nf"
include { check_samplesheet } from '../modules/local/check_samplesheet'

if (params.samplesheet2) { ch_samplesheet = file(params.samplesheet2) } else { exit 1, 'Samplesheet2 not specified!' }


/**
 * Project new data onto the dataset using scANVI
 */
workflow  add_additional_datasets {

    take:
        reference_atlas_h5ad
        reference_scanvi_h5ad
        reference_scanvi_model

    main:

    ch_samples = Channel.from(check_samplesheet(ch_samplesheet.toString()))

    SCQC(
        [
            file("${baseDir}/modules/local/scqc/scqc-notebook.py", checkIfExists: true),
            file("${baseDir}/modules/local/scqc/qc_plots.py", checkIfExists: true)
        ],
        ch_samples
    )
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())

    INTEGRATE_INTO_ATLAS(
        Channel.value([ [id: 'integrate_into_atlas'], file("$baseDir/analyses/36_add_additional_datasets/36_scvi_mapping.py") ]),
        ch_samples.map{meta, adata -> [
            "dataset_path": ".",
            "dataset_id": meta.id,
            "reference_atlas": "full_atlas_annotated.h5ad",
            "reference_scanvi_h5ad": "full_atlas_hvg_integrated_scvi_integrated_scanvi.h5ad",
            "reference_scanvi_model": "full_atlas_hvg_integrated_scvi_scanvi_model",
            "gene_symbol_table": "gene_symbol_dict.csv"
        ]},
        SCQC.out.adata.flatMap{ id, adata -> adata}.mix(
            reference_atlas_h5ad,
            reference_scanvi_h5ad,
            reference_scanvi_model,
            Channel.fromPath("${baseDir}/tables/gene_symbol_dict.csv"),
            ch_samplesheet
        ).collect()
    )

    UPDATE_ANNOTATION(
        [ [id: 'integrate_into_atlas'], file("$baseDir/analyses/36_add_additional_datasets/36b_update_annotation.py") ],
        [
            "dataset_path": ".",
            "reference_atlas": "full_atlas_annotated.h5ad",
            "reference_scanvi_h5ad": "full_atlas_hvg_integrated_scvi_integrated_scanvi.h5ad",
        ],
        INTEGRATE_INTO_ATLAS.out.artifacts.collect().mix(
            reference_atlas_h5ad,
            reference_scanvi_h5ad,
        ).collect()
    )

    // emit:
    //     full_atlas_merged = UPDATE_ANNOTATION.out.artifacts


}

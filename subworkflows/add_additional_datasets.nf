include { SCQC } from "../modules/local/scqc/main"
include { SCQC_MERGE_STATS } from "../modules/local/scqc_merge_stats/main.nf"
include { JUPYTERNOTEBOOK as INTEGRATE_INTO_ATLAS } from "../modules/local/jupyternotebook/main.nf"

/**
 * Project new data onto the dataset using scANVI
 */
workflow  add_additional_datasets {

    take:
        ch_samples
        reference_atlas_h5ad
        reference_scanvi_h5ad
        reference_scanvi_model

    main:

    SCQC(
        [
            file("${baseDir}/modules/local/scqc/scqc-notebook.py", checkIfExists: true),
            file("${baseDir}/modules/local/scqc/qc_plots.py", checkIfExists: true)
        ],
        ch_samples
    )
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())

    INTEGRATE_INTO_ATLAS(
        [ [id: 'integrate_into_atlas'], file("$baseDir/analyses/36_add_additional_datasets/36_scvi_mapping.py") ],
        [
            "dataset_path": ".",
            "samplesheet": "samplesheet_scrnaseq_preprocessing2.csv",
            "reference_atlas": "full_atlas_annotated.h5ad",
            "reference_scanvi_h5ad": "full_atlas_hvg_integrated_scvi_integrated_scanvi.h5ad",
            "reference_scanvi_model": "full_atlas_hvg_integrated_scvi_scanvi_model",
            "gene_symbol_table": "gene_symbol_dict.csv"
        ],
        SCQC.out.adata.flatMap{ id, adata -> adata}.mix(
            Channel.from(reference_atlas_h5ad),
            Channel.from(reference_scanvi_h5ad),
            Channel.from(reference_scanvi_model),
            Channel.fromPath("${baseDir}/tables/gene_symbol_dict.csv"),
            Channel.fromPath("${baseDir}/tables/samplesheet_scrnaseq_preprocessing2.csv")
        ).collect()
    )

    emit:
        full_atlas_merged = INTEGRATE_INTO_ATLAS.out.artifacts


}

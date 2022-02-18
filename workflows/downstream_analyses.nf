#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { check_samplesheet } from '../modules/local/check_samplesheet'
include { add_additional_datasets } from "../subworkflows/add_additional_datasets.nf"
include { de_analysis } from "../subworkflows/de_analysis.nf"
include { scissor } from "../subworkflows/scissor.nf"
include { cell2cell } from "../subworkflows/cell2cell.nf"
include { infercnv } from "../subworkflows/infercnv.nf"
include { plots_and_comparisons } from "../subworkflows/plots_and_comparisons.nf"
include { JUPYTERNOTEBOOK as STRATIFY_PATIENTS } from "../modules/local/jupyternotebook/main.nf"

workflow downstream_analyses {
    assert params.atlas: "Atlas h5ad file not specified!"

    ch_samples = Channel.from(check_samplesheet(params.additional_input, baseDir))
    reference_atlas = file(params.atlas, checkIfExists: true)
    reference_scanvi_h5ad = file(params.reference_scanvi_h5ad, checkIfExists: true)
    reference_scanvi_model = file(params.reference_scanvi_model, checkIfExists: true)

    add_additional_datasets(
        ch_samples,
        reference_atlas,
        reference_scanvi_h5ad,
        reference_scanvi_model
    )

    final_atlas = add_additional_datasets.out.full_atlas_merged

    STRATIFY_PATIENTS(
        Channel.value([
            [id: 'stratify_patients'],
            file("${baseDir}/analyses/38_patient_stratification/38_patient_stratification.py")
        ]),
        final_atlas.map{ it -> ["adata_in": it.name]},
        final_atlas
    )

    de_analysis(final_atlas)
    // scissor(final_atlas)
    cell2cell(final_atlas)
    infercnv(final_atlas)
    plots_and_comparisons(final_atlas, cell2cell.out.adata_cpdb, STRATIFY_PATIENTS.out.artifacts)
}




#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { de_analysis } from "../subworkflows/de_analysis.nf"
include { scissor } from "../subworkflows/scissor.nf"
include { cell2cell } from "../subworkflows/cell2cell.nf"
include { infercnv } from "../subworkflows/infercnv.nf"
include { JUPYTERNOTEBOOK as STRATIFY_PATIENTS } from "../modules/local/jupyternotebook/main.nf"

workflow downstream_analyses {
    assert params.atlas: "Atlas h5ad file not specified!"

    final_atlas = file(params.atlas, checkIfExists: true)

    STRATIFY_PATIENTS(
        Channel.value([
            [id: 'stratify_patients'],
            file("${baseDir}/analyses/38_patient_stratification/38_patient_stratification.py")
        ]),
        ["adata_in": final_atlas.name],
        final_atlas
    )

    de_analysis(final_atlas)
    scissor(final_atlas)
    cell2cell(final_atlas)
    infercnv(final_atlas)
}




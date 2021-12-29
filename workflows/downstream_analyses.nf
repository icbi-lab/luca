#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { de_tumor_normal } from "../subworkflows/de_tumor_normal.nf"
include { scissor } from "../subworkflows/scissor.nf"
include { cell2cell } from "../subworkflows/cell2cell.nf"

workflow downstream_analyses {
    assert params.atlas: "Atlas h5ad file not specified!"

    final_atlas = file(params.atlas, checkIfExists: true)
    de_tumor_normal(final_atlas)
    scissor(final_atlas)
    cell2cell(final_atlas)
}




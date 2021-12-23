#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { de_tumor_normal } from "../subworkflows/de_tumor_normal.nf"
include { scissor } from "../subworkflows/scissor.nf"

workflow downstream_analyses {
    assert params.atlas: "Atlas h5ad file not specified!"

    final_atlas = file(params.atlas)
    de_tumor_normal(final_atlas)
    scissor(final_atlas)
}




#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { build_atlas } from "./workflows/build_atlas.nf"
include { downstream_analyses } from "./workflows/downstream_analyses.nf"

workflow {

    if(params.workflow == "build_atlas") {
        build_atlas()
    } else if (params.workflow == "downstream_analyses") {
        downstream_analyses()
    } else {
        assert False: "Invalid --workflow parameter"
    }

}

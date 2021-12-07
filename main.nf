#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2
assert params.input: "Input samplesheet not specified!"

include { integrate_datasets } from "./workflows/integrate_datasets.nf"
include { annotate_dataset } from "./workflows/annotate_dataset.nf"
include { de_tumor_normal } from "./workflows/de_tumor_normal.nf"
include { scissor } from "./workflows/scissor.nf"

workflow {
    integrate_datasets()
    annotate_dataset(integrate_datasets.out.adata_integrated)
    de_tumor_normal(annotate_dataset.out.final_atlas)
    scissor(integrate_datasets.out.adata_integrated)
}


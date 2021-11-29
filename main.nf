#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2
assert params.input: "Input samplesheet not specified!"

include { integrate_datasets } from "./workflows/integrate_datasets.nf"
include { annotate_dataset } from "./workflows/annotate_dataset.nf"

workflow {
    integrate_datasets()
    annotate_dataset(integrate_datasets.out.adata_integrated)
}


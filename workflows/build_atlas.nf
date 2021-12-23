#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

include { integrate_datasets } from "../subworkflows/integrate_datasets.nf"
include { annotate_dataset } from "../subworkflows/annotate_dataset.nf"

workflow build_atlas {

    assert params.input: "Input samplesheet not specified!"

    integrate_datasets()
    annotate_dataset(integrate_datasets.out.adata_integrated)
}

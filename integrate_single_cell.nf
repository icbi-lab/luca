#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()
assert params.input: "Input samplesheet not specified!"

include { check_samplesheet }  from './modules/local/check_samplesheet' params(params)

include { SCQC } from "./modules/local/scqc/main.nf" addParams(
    options: modules['SCQC']
)
include { SCQC_MERGE_STATS } from "./modules/local/scqc_merge_stats/main.nf" addParams(
    options: modules['SCQC_MERGE_STATS']
)
include { P11_MERGE_ALL } from "./modules/local/analyses/20_integrate_scrnaseq_data.nf" addParams (
    options: modules["P11_MERGE_ALL"]
)


workflow {
    ch_samples = Channel.from(check_samplesheet(params.input))

    SCQC(ch_samples)
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())

    P11_MERGE_ALL(
        Channel.fromPath(params.input),
        SCQC.out.adata.flatMap{ meta, adata -> adata }
    )
}


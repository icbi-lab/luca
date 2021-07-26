#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def modules = params.modules.clone()
assert params.input: "Input samplesheet not specified!"

include { check_samplesheet }  from './modules/local/check_samplesheet' params(params)

include { SCQC } from "./modules/local/scqc/main" addParams(
    options: modules['SCQC']
)
include { SCQC_MERGE_STATS } from "./modules/local/scqc_merge_stats/main.nf" addParams(
    options: modules['SCQC_MERGE_STATS']
)
include { JUPYTERNOTEBOOK as P11_MERGE_ALL } from "./modules/local/jupyternotebook/main.nf" addParams (
    options: modules["P11_MERGE_ALL"]
)
include { SCVI } from "./modules/local/scvi.nf"


workflow {
    ch_samples = Channel.from(check_samplesheet(params.input))

    SCQC(ch_samples)
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())

    P11_MERGE_ALL(
        Channel.value([
            [id: "21_merge_all"],
            file("analyses/20_integrate_scrnaseq_data/21_merge_all.py")
        ]),
        (
            SCQC.out.adata.flatMap{ meta, adata -> adata }.mix(
                Channel.fromPath("tables/samplesheet_scrnaseq_preprocessing.csv")
            )
            .collect()
            .map{ it -> [
                    [
                        samplesheet: "samplesheet_scrnaseq_preprocessing.csv",
                        dataset_path: "."
                    ],
                    it
                ]
            }
        )
    )

    SCVI(
        P11_MERGE_ALL.out.artifacts.collect().map{
            out -> out.findAll{ it -> it.getExtension() == "h5ad" }
        }.view(),
        Channel.from(0, 1)
    )
}


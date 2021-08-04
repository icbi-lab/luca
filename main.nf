#!/usr/bin/env nextflow

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

nextflow.enable.dsl = 2

def modules = params.modules.clone()
def subworkflows = params.subworkflows.clone()
assert params.input: "Input samplesheet not specified!"

include { check_samplesheet }  from './modules/local/check_samplesheet' params(params)

include { SCQC } from "./modules/local/scqc/main" addParams(
    options: modules['SCQC']
)
include { SCQC_MERGE_STATS } from "./modules/local/scqc_merge_stats/main.nf" addParams(
    options: modules['SCQC_MERGE_STATS']
)
include { JUPYTERNOTEBOOK as MERGE_ALL } from "./modules/local/jupyternotebook/main.nf" addParams (
    options: modules["MERGE_ALL"]
)
include { SCVI } from "./modules/local/scvi/main.nf" addParams(
    options: modules["SCVI"]
)
include { SOLO } from "./modules/local/solo/main.nf" addParams(
    options: modules["SOLO"]
)
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_DOUBLET } from "./subworkflows/neighbors_leiden_umap/main.nf" addParams(
    options: subworkflows["NEIGHBORS_LEIDEN_UMAP_DOUBLET"]
)
include { JUPYTERNOTEBOOK as MERGE_SOLO }  from "./modules/local/jupyternotebook/main.nf" addParams (
    options: modules["MERGE_SOLO"]
)
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_COARSE }  from "./modules/local/jupyternotebook/main.nf" addParams (
    options: modules["ANNOTATE_CELL_TYPES_COARSE"]
)
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_CELL_TYPES } from "./subworkflows/neighbors_leiden_umap/main.nf" addParams(
    options: subworkflows["NEIGHBORS_LEIDEN_UMAP_CELL_TYPES"]
)


workflow {
    ch_samples = Channel.from(check_samplesheet(params.input))

    SCQC(ch_samples)
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())

    MERGE_ALL(
        Channel.value([
            [id: "21_merge_all"],
            file("analyses/20_integrate_scrnaseq_data/21_merge_all.py")
        ]),
        [
            samplesheet: "samplesheet_scrnaseq_preprocessing.csv",
            dataset_path: "."
        ],
        SCQC.out.adata.flatMap{ meta, adata -> adata }.mix(
            Channel.fromPath("tables/samplesheet_scrnaseq_preprocessing.csv")
        ).collect()
    )

    SCVI(
        MERGE_ALL.out.artifacts.collect().map{
            out -> out.findAll{ it -> it.getExtension() == "h5ad" }
        },
        Channel.from(0, 1)
    )

    NEIGHBORS_LEIDEN_UMAP_DOUBLET(
        SCVI.out.adata.filter{ it -> !it.baseName.contains("hvg") },
        "X_scVI",
        1.0
    )

    SOLO(
        SCVI.out.adata.filter{ it -> !it.baseName.contains("hvg") },
        SCVI.out.scvi_model.filter{ it -> !it.baseName.contains("hvg") },
        MERGE_ALL.out.artifacts.flatten().filter{
            it -> it.getName() == "obs_all.csv"
        }.splitCsv(header : true).filter{
            it -> it["run_solo"] == "True"
        }.map{ it -> it["sample"] }
    )

    MERGE_SOLO(
        Channel.value([
            [id: "25_merge_solo"],
            file("analyses/20_integrate_scrnaseq_data/25_merge_solo.py")
        ]),
        [
            "adata_path": "integrated_merged_all_all_genes.umap_leiden.h5ad"
        ],
        NEIGHBORS_LEIDEN_UMAP_DOUBLET.out.adata.mix(SOLO.out.doublets).flatten().collect()
    )

    ANNOTATE_CELL_TYPES_COARSE(
        Channel.value([
            [id: "26_annotate_cell_types"],
            file("analyses/20_integrate_scrnaseq_data/26_annotate_cell_types_coarse.py")
        ]),
        [:],
        MERGE_SOLO.out.artifacts
    )

    NEIGHBORS_LEIDEN_UMAP_CELL_TYPES(
        ANNOTATE_CELL_TYPES_COARSE.out.artifacts.flatten(),
        "X_scVI",
        1.0
    )

}


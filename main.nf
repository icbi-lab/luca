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
include { SCVI as SCVI_SEED } from "./modules/local/scvi/main.nf" addParams(
    options: modules["SCVI_SEED"]
)
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_SEED } from "./subworkflows/neighbors_leiden_umap/main.nf" addParams(
    options: subworkflows["NEIGHBORS_LEIDEN_UMAP_SEED"]
)
include { JUPYTERNOTEBOOK as ANNOTATE_SEED } from "./modules/local/jupyternotebook/main.nf" addParams (
    options: modules["ANNOTATE_SEED"]
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



// TODO: Enable "seed annotation" and use SCANVI (SCVI fails to integrate smartseq2 data)
workflow {
    ch_samples = Channel.from(check_samplesheet(params.input, baseDir))

    SCQC(ch_samples)
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())

    // SEED annotation (manually annotate two datasets (One 10x and one Smartseq)
    // in order to use the scANVI algorithm for the integration which has been shown
    // to outperform scVI)
    SCVI_SEED(
        SCQC.out.adata.map{ meta, adata -> [meta.id, adata] }.filter{
            id, adata -> {
                id.equals("Maynard_Bivona_2020_NSCLC") || id.equals("Lambrechts_2018_LUAD_6653")
            }
       },
       1,
       "sample"
    )
    NEIGHBORS_LEIDEN_UMAP_SEED(SCVI_SEED.out.adata, "X_scVI", 1.0)
    ANNOTATE_SEED(
        SCVI_SEED.out.adata.map{ id, adata -> [
            ["id": id],
            file("${baseDir}/analyses/10_seed_annotations/annotate_${id.toLowerCase()}.py")
        ]},
        NEIGHBORS_LEIDEN_UMAP_SEED.out.adata.map{ id, adata -> [
            adata_qc: "${id}.qc.h5ad",
            adata_scvi: adata.name
        ]},
        SCVI_SEED.out.adata.mix(NEIGHBORS_LEIDEN_UMAP_SEED.out.adata).groupTuple().map{
             id, files -> files
        }
    )

    // MERGE and INTEGRATE all datasets



    // MERGE_ALL(
    //     Channel.value([
    //         [id: "21_merge_all"],
    //         file("${baseDir}/analyses/20_integrate_scrnaseq_data/21_merge_all.py")
    //     ]),
    //     [
    //         samplesheet: "samplesheet_scrnaseq_preprocessing.csv",
    //         dataset_path: ".",
    //         gene_symbol_table: "gene_symbol_dict.csv"
    //     ],
    //     SCQC.out.adata.flatMap{ meta, adata -> adata }.mix(
    //         Channel.fromPath("${baseDir}/tables/samplesheet_scrnaseq_preprocessing.csv"),
    //         Channel.fromPath("${baseDir}/tables/gene_symbol_dict.csv")
    //     ).collect()
    // )

    // SCVI(
    //     MERGE_ALL.out.artifacts.collect().map{
    //         out -> out.findAll{ it -> it.getExtension() == "h5ad" }
    //     },
    //     Channel.from(0, 1)
    // )

    // NEIGHBORS_LEIDEN_UMAP_DOUBLET(
    //     SCVI.out.adata.filter{ it -> it.baseName.contains("hvg") },
    //     "X_scVI",
    //     1.0
    // )

    // SOLO(
    //     SCVI.out.adata.filter{ it -> it.baseName.contains("hvg") },
    //     SCVI.out.scvi_model.filter{ it -> it.baseName.contains("hvg") },
    //     MERGE_ALL.out.artifacts.flatten().filter{
    //         it -> it.getName() == "obs_all.csv"
    //     }.splitCsv(header : true).filter{
    //         it -> it["run_solo"] == "True"
    //     }.map{ it -> it["sample"] }
    // )

    // MERGE_SOLO(
    //     Channel.value([
    //         [id: "25_merge_solo"],
    //         file("${baseDir}/analyses/20_integrate_scrnaseq_data/25_merge_solo.py")
    //     ]),
    //     [
    //         "adata_path": "integrated_merged_all_hvg.umap_leiden.h5ad"
    //     ],
    //     NEIGHBORS_LEIDEN_UMAP_DOUBLET.out.adata.mix(SOLO.out.doublets).flatten().collect()
    // )

    // ANNOTATE_CELL_TYPES_COARSE(
    //     Channel.value([
    //         [id: "26_annotate_cell_types"],
    //         file("${baseDir}/analyses/20_integrate_scrnaseq_data/26_annotate_cell_types_coarse.py")
    //     ]),
    //     [:],
    //     MERGE_SOLO.out.artifacts
    // )

    // NEIGHBORS_LEIDEN_UMAP_CELL_TYPES(
    //     ANNOTATE_CELL_TYPES_COARSE.out.artifacts.flatten().filter( it -> !it.baseName.contains("cell_type_coarse")),
    //     "X_scVI",
    //     Channel.from(0.5, 0.75, 1.0, 1.5)
    // )

}


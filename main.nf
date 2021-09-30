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
include { SCANVI } from "./modules/local/scvi/main.nf" addParams(
    options: modules["SCANVI"]
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
include { JUPYTERNOTEBOOK as ANNOTATE_CELL_TYPES_FINE }  from "./modules/local/jupyternotebook/main.nf" addParams (
    options: modules["ANNOTATE_CELL_TYPES_FINE"]
)
include { JUPYTERNOTEBOOK as PREPARE_CELLXGENE }  from "./modules/local/jupyternotebook/main.nf" addParams (
    options: modules["PREPARE_CELLXGENE"]
)
include { H5AD_TO_SEURAT }  from "./modules/local/scconversion/main.nf" addParams(
    options: ["publish_dir": "test-conversion"]
)
include { SEURAT_TO_SCE }  from "./modules/local/scconversion/main.nf" addParams(
    options: ["publish_dir": "test-conversion"]
)



// TODO: Enable "seed annotation" and use SCANVI (SCVI fails to integrate smartseq2 data)
workflow {
    ch_samples = Channel.from(check_samplesheet(params.input, baseDir))

    SCQC(ch_samples)
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())

    // SEED annotation (manually annotate two datasets (One 10x and one Smartseq)
    // in order to use the scANVI algorithm for the integration which has been shown
    // to outperform scVI)
    ch_seed_ids = Channel.from("Maynard_Bivona_2020_NSCLC", "Lambrechts_2018_LUAD_6653")
    SCVI_SEED(
       ch_seed_ids.join(SCQC.out.adata),
       1,
       ["sample", "sample", null]
    )
    NEIGHBORS_LEIDEN_UMAP_SEED(SCVI_SEED.out.adata, "X_scVI", 1.0)
    ch_seed_scvi = SCQC.out.adata.join(NEIGHBORS_LEIDEN_UMAP_SEED.out.adata)
    ANNOTATE_SEED(
        ch_seed_scvi.map{id, adata1, adata2 -> [
            ["id": id],
            file("${baseDir}/analyses/10_seed_annotations/annotate_${id.toLowerCase()}.py")
        ]},
        ch_seed_scvi.map{ id, adata_qc, adata_scvi -> [
            adata_qc: adata_qc.name,
            adata_scvi: adata_scvi.name
        ]},
        ch_seed_scvi.map{
             id, adata_qc, adata_scvi -> [adata_qc, adata_scvi]
        }
    )

    // MERGE and INTEGRATE all datasets
    MERGE_ALL(
        channel.value([
            [id: "21_merge_all"],
            file("${baseDir}/analyses/20_integrate_scrnaseq_data/21_merge_all.py")
        ]),
        [
            samplesheet: "samplesheet_scrnaseq_preprocessing.csv",
            dataset_path: ".",
            dataset_path_annotated: ".",
            gene_symbol_table: "gene_symbol_dict.csv"
        ],
        SCQC.out.adata.flatMap{ id, adata -> adata }.mix(
            ANNOTATE_SEED.out.artifacts
        ).mix(
            Channel.fromPath("${baseDir}/tables/samplesheet_scrnaseq_preprocessing.csv"),
            Channel.fromPath("${baseDir}/tables/gene_symbol_dict.csv")
        ).collect()
    )

    SCVI(
        MERGE_ALL.out.artifacts.collect().map{
            out -> ["all", out.findAll{ it -> it.getExtension() == "h5ad" }]
        },
        Channel.from(0, 1),
        ["batch", "dataset", null]
    )

    SCANVI(
        SCVI.out.adata.join(SCVI.out.scvi_model),
        "batch",
        "cell_type"
    )

    // use HVG version for downstream analysis. We just keep the version with
    // all genes in case we need to run DE analysis with scVI.
    ch_scvi_hvg = SCVI.out.adata.filter{ id, adata -> adata.baseName.contains("hvg") }
    ch_scvi_hvg_model = SCVI.out.scvi_model.filter{ id, adata -> adata.baseName.contains("hvg") }
    ch_scanvi_hvg = SCANVI.out.adata.filter{ id, adata -> adata.baseName.contains("hvg") }
    ch_scanvi_hvg_model = SCANVI.out.scvi_model.filter{ id, adata -> adata.baseName.contains("hvg") }

    NEIGHBORS_LEIDEN_UMAP_DOUBLET(
        ch_scanvi_hvg,
        "X_scANVI",
        1.0
    )

    SOLO(
        ch_scvi_hvg,
        ch_scvi_hvg_model,
        MERGE_ALL.out.artifacts.flatten().filter{
            it -> it.getName() == "obs_all.csv"
        }.splitCsv(header : true).filter{
            it -> it["run_solo"] == "True"
        }.map{ it -> it["sample"] }
    )

    MERGE_SOLO(
        Channel.value([
            [id: "25_merge_solo"],
            file("${baseDir}/analyses/20_integrate_scrnaseq_data/25_merge_solo.py")
        ]),
        [
            "adata_path": "all.umap_leiden.h5ad"
        ],
        NEIGHBORS_LEIDEN_UMAP_DOUBLET.out.adata.map{ id, adata -> adata}.mix(
            SOLO.out.doublets
        ).flatten().collect()
    )

    ANNOTATE_CELL_TYPES_COARSE(
        Channel.value([
            [id: "27_annotate_cell_types"],
            file("${baseDir}/analyses/20_integrate_scrnaseq_data/27_annotate_cell_types_coarse.py")
        ]),
        [:],
        MERGE_SOLO.out.artifacts
    )

    ch_adata_annotated_by_cell_type = ANNOTATE_CELL_TYPES_COARSE.out.artifacts.flatten().filter(
        it -> !it.baseName.contains("cell_type_coarse")
    )
    ch_adata_annotated = ANNOTATE_CELL_TYPES_COARSE.out.artifacts.flatten().filter(
        it -> it.baseName.contains("cell_type_coarse")
    )

    NEIGHBORS_LEIDEN_UMAP_CELL_TYPES(
        ch_adata_annotated_by_cell_type.map{ adata -> [adata.baseName, adata] },
        "X_scANVI",
        Channel.from(0.5, 0.75, 1.0, 1.5)
    )

    ANNOTATE_CELL_TYPES_FINE(
        Channel.value([
            [id: "29_annotate_cell_types_fine"],
            file("${baseDir}/analyses/20_integrate_scrnaseq_data/29_annotate_cell_types_fine.py")
        ]),
        [
            "input_dir": '.',
            "main_adata": 'adata_cell_type_coarse.h5ad'
        ],
        NEIGHBORS_LEIDEN_UMAP_CELL_TYPES.out.adata.map{ id, adata -> adata }.mix(
            ch_adata_annotated
        ).collect()
    )

    PREPARE_CELLXGENE(
        Channel.value([
            [id: "zz_prepare_cellxgene"],
            file("${baseDir}/analyses/zz_cellxgene/stats_and_cellxgene.py")
        ]),
        ["adata_in": "adata_annotated_fine.h5ad"],
        ANNOTATE_CELL_TYPES_FINE.out.artifacts
    )

    H5AD_TO_SEURAT(
        Channel.value(['organoids', file('/data/projects/2017/Organoids-ICBI/zenodo/scrnaseq/03_scvi/adata_integrated.h5ad')])
    )
    SEURAT_TO_SCE(H5AD_TO_SEURAT.out.h5seurat)

}


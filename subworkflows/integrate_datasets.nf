include { check_samplesheet } from '../modules/local/check_samplesheet'

include { SCQC } from "../modules/local/scqc/main"
include { SCQC_MERGE_STATS } from "../modules/local/scqc_merge_stats/main.nf"

include { JUPYTERNOTEBOOK as MERGE_ALL } from "../modules/local/jupyternotebook/main.nf"
include { SCVI } from "../modules/local/scvi/main.nf"
include { SCANVI } from "../modules/local/scvi/main.nf"
include { SOLO } from "../modules/local/solo/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_DOUBLET } from "./neighbors_leiden_umap.nf"
include { JUPYTERNOTEBOOK as MERGE_SOLO }  from "../modules/local/jupyternotebook/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_NODOUBLET } from "./neighbors_leiden_umap.nf"

if (params.samplesheet) { ch_samplesheet = file(params.samplesheet) } else { exit 1, 'Samplesheet not specified!' }


/**
 * Integrate individual datasets into a single-cell atlas
 *   - Perfom QC on individual datasets
 *   - Perform manual "seed" annotation of two datasets (one SS2, one 10x)
 *   - Perform data integration using SCANVI
 *   - Call doublets using solo
 */
workflow integrate_datasets {


    main:

    ch_samples = Channel.from(check_samplesheet(ch_samplesheet.toString()))

    SCQC(
        [
            file("${baseDir}/modules/local/scqc/scqc-notebook.py", checkIfExists: true),
            file("${baseDir}/modules/local/scqc/qc_plots.py", checkIfExists: true)
        ],
        ch_samples
    )
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())


    // MERGE and INTEGRATE all datasets
    MERGE_ALL(
        Channel.value([
            [id: "21_merge_all"],
            file("${baseDir}/analyses/20_integrate_scrnaseq_data/21_merge_all.py")
        ]),
        [
            samplesheet: ch_samplesheet.toString(),
            gene_symbol_table: "gene_symbol_dict.csv"
        ],
        SCQC.out.adata.flatMap{ id, adata -> adata }.mix(
            Channel.value(ch_samplesheet),
            Channel.fromPath("${baseDir}/tables/gene_symbol_dict.csv")
        ).collect()
    )
    

    ch_adata_merged = MERGE_ALL.out.artifacts.collect().map{
        out -> ["all", out.findAll{ it -> it.getExtension() == "h5ad" }]
    }

    SCVI(
        ch_adata_merged,
        1, // 1 = use HVG
        ["batch", "dataset", null]
    )

    /*
    SCANVI(
        SCVI.out.adata.join(SCVI.out.scvi_model),
        "batch",
        "cell_type"
    )
    ch_scvi_hvg = SCVI.out.adata
    ch_scvi_hvg_model = SCVI.out.scvi_model
    ch_scanvi_hvg = SCANVI.out.adata
    ch_scanvi_hvg_model = SCANVI.out.scvi_model

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
            "adata_path": "all.umap_leiden.h5ad",
            // this is to re-integrate all genes (not only HVG)
            "adata_merged": "merged_all.h5ad"
        ],
        NEIGHBORS_LEIDEN_UMAP_DOUBLET.out.adata.map{ id, adata -> adata}.mix(
            SOLO.out.doublets
        ).mix(ch_adata_merged.map{ id, it -> it}).flatten().collect()
    )

    //re-compute neighbors, leiden, umap after doublet filtering.
    ch_adata_doublet_filtered = MERGE_SOLO.out.artifacts.filter{
        it -> it.baseName.contains("doublet_filtered")
    }.map{ it -> [it.baseName, it] }
    NEIGHBORS_LEIDEN_UMAP_NODOUBLET(
        ch_adata_doublet_filtered,
        "X_scANVI",
        1.0
    )

    emit:
        adata_integrated = NEIGHBORS_LEIDEN_UMAP_NODOUBLET.out.adata.map{ meta, ad -> ad }
    */
}

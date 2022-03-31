include { check_samplesheet } from '../modules/local/check_samplesheet'

include { SCQC } from "../modules/local/scqc/main"
include { SCQC_MERGE_STATS } from "../modules/local/scqc_merge_stats/main.nf"
include { SCVI as SCVI_SEED } from "../modules/local/scvi/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_SEED } from "./neighbors_leiden_umap.nf"
include { JUPYTERNOTEBOOK as ANNOTATE_SEED } from "../modules/local/jupyternotebook/main.nf"


include { JUPYTERNOTEBOOK as MERGE_ALL } from "../modules/local/jupyternotebook/main.nf"
include { SCVI } from "../modules/local/scvi/main.nf"
include { SCANVI } from "../modules/local/scvi/main.nf"
include { SOLO } from "../modules/local/solo/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_DOUBLET } from "./neighbors_leiden_umap.nf"
include { JUPYTERNOTEBOOK as MERGE_SOLO }  from "../modules/local/jupyternotebook/main.nf"
include { NEIGHBORS_LEIDEN_UMAP as NEIGHBORS_LEIDEN_UMAP_NODOUBLET } from "./neighbors_leiden_umap.nf"


/**
 * Integrate individual datasets into a single-cell atlas
 *   - Perfom QC on individual datasets
 *   - Perform manual "seed" annotation of two datasets (one SS2, one 10x)
 *   - Perform data integration using SCANVI
 *   - Call doublets using solo
 */
workflow integrate_datasets {


    main:

    ch_samples = Channel.from(check_samplesheet("${baseDir}/tables/samplesheet_scrnaseq_preprocessing.csv", baseDir))
    SCQC(
        [
            file("${baseDir}/modules/local/scqc/scqc-notebook.py", checkIfExists: true),
            file("${baseDir}/modules/local/scqc/qc_plots.py", checkIfExists: true)
        ],
        ch_samples
    )
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

    ch_adata_merged = MERGE_ALL.out.artifacts.collect().map{
        out -> ["all", out.findAll{ it -> it.getExtension() == "h5ad" }]
    }

    SCVI(
        ch_adata_merged,
        1, // 1 = use HVG
        ["batch", "dataset", null]
    )
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
}

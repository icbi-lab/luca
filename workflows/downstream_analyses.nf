#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { de_analysis } from "../subworkflows/de_analysis.nf"
include { scissor } from "../subworkflows/scissor.nf"
include { infercnv } from "../subworkflows/infercnv.nf"
include { plots_and_comparisons } from "../subworkflows/plots_and_comparisons.nf"
include { JUPYTERNOTEBOOK as STRATIFY_PATIENTS } from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as NEUTROPHIL_SUBCLUSTERING } from "../modules/local/jupyternotebook/main.nf"

workflow downstream_analyses {
    assert params.atlas: "Atlas h5ad file not specified!"

    final_atlas = Channel.fromPath(params.atlas)

    NEUTROPHIL_SUBCLUSTERING(
        Channel.value([
            [id: 'neutrophil_subclustering'],
            file("${baseDir}/analyses/37_subclustering/37_neutrophil_subclustering.py")
        ]),
        final_atlas.map{ it -> ["adata_in": it.name, "neutro_clustering": "neutrophil_clustering.csv"]},
        final_atlas.mix(Channel.fromPath("${baseDir}/tables/neutrophil_clustering.csv")).collect()
    )

    atlas_neutro_clusters = NEUTROPHIL_SUBCLUSTERING.out.artifacts.flatten().filter{ it -> it.baseName.equals("full_atlas_neutrophil_clusters") }
    neutro_clusters = NEUTROPHIL_SUBCLUSTERING.out.artifacts.flatten().filter{ it -> it.baseName.equals("adata_neutrophil_clusters") }

    STRATIFY_PATIENTS(
        Channel.value([
            [id: 'stratify_patients'],
            file("${baseDir}/analyses/38_patient_stratification/38_patient_stratification.py")
        ]),
        final_atlas.map{ it -> ["adata_in": it.name]},
        final_atlas
    )

    de_analysis(final_atlas, STRATIFY_PATIENTS.out.artifacts)

    // scissor(final_atlas)
    infercnv(final_atlas)
    plots_and_comparisons(final_atlas, STRATIFY_PATIENTS.out.artifacts)
}




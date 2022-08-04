#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { de_analysis } from "../subworkflows/de_analysis.nf"
include { scissor } from "../subworkflows/scissor.nf"
include { infercnv } from "../subworkflows/infercnv.nf"
include { plots_and_comparisons } from "../subworkflows/plots_and_comparisons.nf"
include { JUPYTERNOTEBOOK as STRATIFY_PATIENTS } from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as NEUTROPHIL_SUBCLUSTERING } from "../modules/local/jupyternotebook/main.nf"

workflow downstream_analyses {
    assert params.build_atlas_dir: "Atlas h5ad file not specified!"

    // Get input data from upstream `build_atlas` workflow
    core_atlas = Channel.fromPath(
        "${params.build_atlas_dir}/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad",
         checkIfExists: true
    )
    core_atlas_epithelial_cells = Channel.fromPath(
        "${params.build_atlas_dir}/annotate_datasets/35_final_atlas/artifacts/epithelial_cells_annotated.h5ad",
         checkIfExists: true
    )
    core_atlas_tumor_cells = Channel.fromPath(
        "${params.build_atlas_dir}/annotate_datasets/33_cell_types_epi/artifacts/adata_tumor.h5ad",
         checkIfExists: true
    )
    extended_atlas = Channel.fromPath(
        "${params.build_atlas_dir}/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad",
         checkIfExists: true
    )


    NEUTROPHIL_SUBCLUSTERING(
        Channel.value([
            [id: 'neutrophil_subclustering'],
            file("${baseDir}/analyses/37_subclustering/37_neutrophil_subclustering.py")
        ]),
        extended_atlas.map{ it -> ["adata_in": it.name]},
        extended_atlas
    )
    atlas_neutro_clusters = NEUTROPHIL_SUBCLUSTERING.out.artifacts.flatten().filter{ it -> it.baseName.equals("full_atlas_neutrophil_clusters") }
    neutro_clusters = NEUTROPHIL_SUBCLUSTERING.out.artifacts.flatten().filter{ it -> it.baseName.equals("adata_neutrophil_clusters") }

    STRATIFY_PATIENTS(
        Channel.value([
            [id: 'stratify_patients'],
            file("${baseDir}/analyses/38_patient_stratification/38_patient_stratification.py")
        ]),
        extended_atlas.map{ it -> ["adata_in": it.name]},
        extended_atlas
    )
    patient_stratification_table = STRATIFY_PATIENTS.out.artifacts.flatten().filter{ it -> it.name.equals("patient_stratification.csv") }
    patient_stratification_adata_immune = STRATIFY_PATIENTS.out.artifacts.flatten().filter{ it -> it.name.equals("adata_immune.h5ad") }
    patient_stratification_adata_tumor_subtypes = STRATIFY_PATIENTS.out.artifacts.flatten().filter{ it -> it.name.equals("adata_tumor_subtypes.h5ad") }

    de_analysis(core_atlas, extended_atlas, patient_stratification_table)
    de_result_tumor_cells = de_analysis.out.immune_infiltration.mix(
        de_analysis.out.luad_lusc
    ).flatten().filter{ it -> it.baseName.contains("tumor_cells") }

    // scissor(extended_atlas)
    // infercnv(extended_atlas, patient_stratification_table)
    plots_and_comparisons(
        extended_atlas,
        neutro_clusters,
        core_atlas,
        core_atlas_epithelial_cells,
        core_atlas_tumor_cells,
        patient_stratification_table,
        patient_stratification_adata_immune,
        patient_stratification_adata_tumor_subtypes,
        de_result_tumor_cells
    )
}




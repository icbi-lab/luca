#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { de_analysis } from "../subworkflows/de_analysis.nf"
include { scissor } from "../subworkflows/scissor.nf"
include { plots_and_comparisons } from "../subworkflows/plots_and_comparisons.nf"
include { JUPYTERNOTEBOOK as STRATIFY_PATIENTS } from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as STRATIFY_PATIENTS_SAMPLING_LOCATION } from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as NEUTROPHIL_SUBCLUSTERING } from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as EXPORT_ATLAS } from "../modules/local/jupyternotebook/main.nf"

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

    STRATIFY_PATIENTS_SAMPLING_LOCATION(
        Channel.value([
            [id: "stratify_patients_sampling_location"],
            file("${baseDir}/analyses/38_patient_stratification/38b_sampling_location_lambrechts.py")
        ]),
        extended_atlas.map{ it -> ["adata_in": it.name]},
        extended_atlas
    )

    de_analysis(core_atlas, extended_atlas, neutro_clusters, patient_stratification_table)
    de_result_tumor_cells = de_analysis.out.immune_infiltration.mix(
        de_analysis.out.luad_lusc
    ).flatten().filter{ it -> it.baseName.contains("tumor_cells") }
    de_result_tan_nan = de_analysis.out.tan_nan.flatten().filter{ it -> it.name ==~ /tan_nan.*_DESeq2_result.tsv/ }
    de_result_neutro_clusters = de_analysis.out.neutro_clusters.flatten().filter{ it -> it.name ==~ /neutrophil_subclusters.*_DESeq2_result.tsv/ }

    scissor(extended_atlas)
    plots_and_comparisons(
        extended_atlas,
        neutro_clusters,
        core_atlas,
        core_atlas_epithelial_cells,
        core_atlas_tumor_cells,
        patient_stratification_table,
        patient_stratification_adata_immune,
        patient_stratification_adata_tumor_subtypes,
        de_result_tumor_cells,
        de_result_tan_nan,
        de_result_neutro_clusters
    )

    ch_symbol_to_ensembl = Channel.fromPath("${baseDir}/tables/symbol_to_ensembl.csv")
    ch_export_atlas_extended = extended_atlas.concat(atlas_neutro_clusters, ch_symbol_to_ensembl).collect()
    ch_export_atlas_core = core_atlas.concat(ch_symbol_to_ensembl).collect()
    ch_export_atlas_extended_params = ch_export_atlas_extended.map{
        extended_atlas, neutro_atlas, symbol_to_ensembl -> [
            "id": "extended_atlas",
            "atlas": extended_atlas.name,
            "neutrophil_atlas": neutro_atlas.name,
            "title": "The single-cell lung cancer atlas (LuCA) -- extended atlas",
            "output_filename": "extended_atlas_cellxgene_schema.h5ad",
            "symbol_to_ensembl": symbol_to_ensembl.name
        ]
    }
    ch_export_atlas_core_params = ch_export_atlas_core.map{
        core_atlas, symbol_to_ensembl -> [
            "id": "core_atlas",
            "atlas": core_atlas.name,
            "neutrophil_atlas": "None",
            "title": "The single-cell lung cancer atlas (LuCA) -- core atlas",
            "output_filename": "core_atlas_cellxgene_schema.h5ad",
            "symbol_to_ensembl": symbol_to_ensembl.name
        ]
    }
    EXPORT_ATLAS(
        Channel.value([
            [id: "export_atlas"],
            file("${baseDir}/analyses/80_export_atlas/81_export_atlas.py")
        ]),
        ch_export_atlas_extended_params.concat(ch_export_atlas_core_params),
        ch_export_atlas_extended.concat(ch_export_atlas_core)
    )
}




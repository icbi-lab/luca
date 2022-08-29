include { JUPYTERNOTEBOOK as PREPARE_FOR_DE } from "../modules/local/jupyternotebook/main.nf"
include { deseq2_analysis as de_analysis_tumor_normal;
          deseq2_analysis as de_analysis_luad_lusc;
          deseq2_analysis as de_analysis_early_advanced;
          deseq2_analysis as de_analysis_immune_infiltration;
          deseq2_analysis as de_analysis_tumor_cell_types;
          deseq2_analysis as de_analysis_tan_nan;
          deseq2_analysis as de_analysis_neutro_clusters;
        } from "../modules/local/scde/main.nf"
include { SPLIT_ANNDATA } from "../modules/local/scconversion/main.nf"


 workflow de_analysis {
    take:
    core_atlas
    final_atlas
    neutro_clusters
    patient_stratification

    main:
    ch_prepare_de_input = final_atlas.concat(patient_stratification).collect()
    PREPARE_FOR_DE(
        Channel.value([
            [id: "prepare_for_de"],
            file("${baseDir}/analyses/40_de_analysis/41_prepare_de_analysis.py")
        ]),
        ch_prepare_de_input.map{final_atlas, patient_table -> [
            "input_adata": final_atlas.name,
            "patient_stratification": patient_table.name
        ]},
        ch_prepare_de_input
    )
    ch_prepare_adata = PREPARE_FOR_DE.out.artifacts.flatten().map { it -> [it.baseName, it] }
    ch_adata_tumor_normal = ch_prepare_adata.filter{ id, adata -> id == "adata_tumor_normal"}
    ch_adata_primary_tumor = ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"}
    ch_adata_neutro_clusters = neutro_clusters.map{ it -> [it.baseName, it]}

    de_analysis_tumor_normal(
        "tumor_normal",
        ch_adata_tumor_normal,
        "origin",
        ["tumor_primary", "normal_adjacent"],
        "cell_type_major",
        "patient",
        [10, false],
        6, // keep only cell-types with at least 3 paired samples
        "+ patient",
        null
    )
    de_analysis_luad_lusc(
        "luad_lusc_primary_tumor",
        ch_adata_primary_tumor,
        "condition",
        ["LUAD", "LUSC"],
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10 samples
        "+ dataset",
        null
    )
    de_analysis_early_advanced(
        "early_advanced_primary_tumor",
        ch_adata_primary_tumor,
        "tumor_stage",
        ["early", "advanced"],
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10  samples
        "+ dataset",
        null
    )
    de_analysis_immune_infiltration(
        "immune_infiltration_primary_tumor",
        ch_adata_primary_tumor,
        "immune_infiltration",
        "sum2zero",
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10 samples
        "+ dataset",
        null
    )

    /**
     * Find empirical marker genes for tumor cell clusters
     */
    de_analysis_tumor_cell_types(
        "core_atlas_tumor_cell_types",
        core_atlas.map{ it -> ["core_atlas", it]},
        "cell_type_tumor",
        "sum2zero",
        "cell_type",
        "patient",
        [10, true],
        10,
        "+ dataset",
        ["tumor_cells"]
    )

    de_analysis_tan_nan(
        "tan_nan_markers",
        ch_adata_neutro_clusters,
        "cell_type_tan_nan",
        ["TAN", "NAN"],
        "cell_type_coarse",
        "patient",
        [10, true],
        10,
        "+ dataset",
        ["neutrophils"]
    )

    de_analysis_neutro_clusters(
        "neutrophil_subclusters",
        ch_adata_neutro_clusters,
        "cell_type",
        "sum2zero",
        "cell_type_coarse",
        "patient",
        [10, true],
        10,
        "+ dataset",
        ["neutrophils"]
    )

    emit:
    luad_lusc = de_analysis_luad_lusc.out.deseq2_result
    immune_infiltration = de_analysis_immune_infiltration.out.deseq2_result
    tan_nan = de_analysis_tan_nan.out.deseq2_result
    neutro_clusters = de_analysis_neutro_clusters.out.deseq2_result

 }

include { JUPYTERNOTEBOOK as PREPARE_FOR_DE } from "../modules/local/jupyternotebook/main.nf"
include { deseq2_analysis as de_analysis_tumor_normal;
          deseq2_analysis as de_analysis_luad_lusc;
          deseq2_analysis as de_analysis_early_advanced;
          deseq2_analysis as de_analysis_t_desert;
          deseq2_analysis as de_analysis_m_desert;
          deseq2_analysis as de_analysis_b_desert } from "../modules/local/scde/main.nf"


 workflow de_analysis {
    take:
    final_atlas
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

    de_analysis_tumor_normal(
        "tumor_normal",
        ch_prepare_adata.filter{ id, adata -> id == "adata_tumor_normal"},
        "origin",
        [["tumor_primary"], "rest"],
        "cell_type_major",
        "patient",
        [10, false],
        6, // keep only cell-types with at least 3 paired samples
        "+ patient"
    )
    de_analysis_luad_lusc(
        "luad_lusc_primary_tumor",
        ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"},
        "condition",
        [["LUAD"], ["LUSC"]],
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10 samples
        "+ dataset"
    )
    de_analysis_early_advanced(
        "early_advanced_primary_tumor",
        ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"},
        "tumor_stage",
        [["early"], ["advanced"]],
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10  samples
        "+ dataset"
    )
    de_analysis_b_desert(
        "b_desert_primary_tumor",
        ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"},
        "immune_infiltration",
        [["B"], ["desert"]],
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10 samples
        "+ dataset"
    )
    de_analysis_t_desert(
        "t_desert_primary_tumor",
        ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"},
        "immune_infiltration",
        [["T"], ["desert"]],
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10 samples
        "+ dataset"
    )
    de_analysis_m_desert(
        "m_desert_primary_tumor",
        ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"},
        "immune_infiltration",
        [["M"], ["desert"]],
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10 samples
        "+ dataset"
    )

    emit:
    luad_lusc = de_analysis_luad_lusc.out.deseq2_result
    b_desert = de_analysis_b_desert.out.deseq2_result
    t_desert = de_analysis_t_desert.out.deseq2_result
    m_desert = de_analysis_m_desert.out.deseq2_result

 }

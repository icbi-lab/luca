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
    PREPARE_FOR_DE(
        Channel.value([
            [id: "prepare_for_de"],
            file("${baseDir}/analyses/40_de_analysis/41_prepare_de_analysis.py")
        ]),
        final_atlas.map{it -> [
            "input_adata": it.name,
            "patient_stratification": "patient_stratification.csv"
        ]},
        final_atlas.mix(patient_stratification).collect()
    )
    ch_prepare_adata = PREPARE_FOR_DE.out.artifacts.flatten().map { it -> [it.baseName, it] }

    de_analysis_tumor_normal(
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
        ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"},
        "immune_infiltration",
        [["M"], ["desert"]],
        "cell_type_major",
        "patient",
        [10, true],
        10, // keep only cell-types with at least 10 samples
        "+ dataset"
    )

 }

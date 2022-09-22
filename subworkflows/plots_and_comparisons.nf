include {
    JUPYTERNOTEBOOK as NEUTROPHIL_ANALYSIS;
    JUPYTERNOTEBOOK as NEUTROPHIL_ANALYSIS_VELOCYTO;
    JUPYTERNOTEBOOK as STRATIFY_PATIENTS_FIGURES;
    JUPYTERNOTEBOOK as COMPARE_GROUPS;
    JUPYTERNOTEBOOK as COMPARE_GROUPS_PLOTS;
    JUPYTERNOTEBOOK as CELL_TYPE_MARKERS_CORE_ATLAS;
    JUPYTERNOTEBOOK as OVERVIEW_PLOTS_CORE_ATLAS;
    JUPYTERNOTEBOOK as OVERVIEW_PLOTS_EXTENDED_ATLAS;
    JUPYTERNOTEBOOK as CELL_TYPE_SIGNATURES;
    JUPYTERNOTEBOOK as COMPARE_PLATFORMS;
    JUPYTERNOTEBOOK as CPDB_ANALYSIS;
    JUPYTERNOTEBOOK as SCCODA_CONDITION;
} from "../modules/local/jupyternotebook/main.nf"
include {
    RMARKDOWNNOTEBOOK as RESPONSE_TO_ICI;
    RMARKDOWNNOTEBOOK as CPDB_CIRCOS_PLOT;
} from "../modules/local/rmarkdownnotebook/main.nf"

workflow plots_and_comparisons {
    take:
    extended_atlas
    adata_neutrophil_clusters
    core_atlas
    core_atlas_epithelial_cells
    core_atlas_tumor_cells
    patient_stratification
    patient_stratification_adata_immune
    patient_stratification_adata_tumor_subtypes
    deseq2_results
    deseq2_results_tan_nan
    deseq2_results_neutro_clusters

    main:

    ch_neutrophil_analysis_input_files = extended_atlas.concat(
        adata_neutrophil_clusters,
        patient_stratification,
        Channel.fromPath("${baseDir}/tables/gene_annotations/neutro_phenotype_genesets.xlsx"),
        Channel.fromPath("${baseDir}/tables/gene_annotations/neutro_recruitment_chemokines.xlsx"),
        deseq2_results_tan_nan
    ).collect()
    NEUTROPHIL_ANALYSIS(
        Channel.value(
            [[id: 'neutrophil_analysis'], file("${baseDir}/analyses/90_plots_and_comparisons/96a_neutrophils.py")]
        ),
        ch_neutrophil_analysis_input_files.map { adata, adata_n, patient_strat, genesets1, genesets2, deseq2_tan_nan -> [
            "adata_n_path": adata_n.name,
            "adata_path": adata.name,
            "patient_stratification_path": patient_strat.name,
            "neutro_geneset_path": genesets1.name,
            "neutro_recruitment_geneset_path": genesets2.name,
            "deseq_tan_nan_path": deseq2_tan_nan.name
        ]},
        ch_neutrophil_analysis_input_files
    )
    ch_neutro_sigs = NEUTROPHIL_ANALYSIS.out.artifacts.flatten().filter{ it -> it.name.equals("neutro_sigs.csv") }

    ch_velocyto_input_files = adata_neutrophil_clusters.concat(
        Channel.fromPath("${baseDir}/data/11_own_datasets/velocyto/")
    ).collect()
    NEUTROPHIL_ANALYSIS_VELOCYTO(
        Channel.value(
            [[id: 'neutrophil_analysis_velocyto'], file("${baseDir}/analyses/90_plots_and_comparisons/96b_neutrophils_velocyto.py")]
        ),
        ch_velocyto_input_files.map { adata_n, velocyto_dir -> [
            "adata_n_path": adata_n.name,
            "velocyto_dir": velocyto_dir.name
        ]},
        ch_velocyto_input_files
    )


    ch_patient_stratification_figures_input_files = patient_stratification.concat(
        patient_stratification_adata_immune,
        patient_stratification_adata_tumor_subtypes
    ).collect()
    STRATIFY_PATIENTS_FIGURES(
        Channel.value(
            [[id: 'patient_stratification_figures'], file("${baseDir}/analyses/90_plots_and_comparisons/93_patient_stratification_figures.py")]
        ),
        ch_patient_stratification_figures_input_files.map { pat_table, ad_immune, ad_tumor_subtypes -> [
            "patient_stratification_path": pat_table.name,
            "ad_immune_path": ad_immune.name,
            "ad_tumor_subtypes_path": ad_tumor_subtypes.name
        ]},
        ch_patient_stratification_figures_input_files
    )


    COMPARE_GROUPS(
        Channel.value(
            [[id: 'compare_groups'], file("${baseDir}/analyses/90_plots_and_comparisons/91_compare_groups.py")]
        ),
        Channel.from([
            // "tumor_normal",
            "patient_infiltration_status",
            "patient_immune_infiltration",
            "patient_immune_infiltration_treatment_coding",
            "patient_immune_infiltration_treatment_coding_random",
            "luad_lusc",
            "early_advanced",
        ]).map {
            it -> [
                "comparison": it,
                "adata_in": "full_atlas_merged.h5ad",
                "stratification_csv": "patient_stratification.csv"
            ]
        },
        extended_atlas.mix(patient_stratification).collect()
    )

    COMPARE_GROUPS_PLOTS(
        Channel.value(
            [[id: 'compare_groups_plots'], file("${baseDir}/analyses/90_plots_and_comparisons/92_compare_groups_plots.py")]
        ),
        [
            "path_prefix": "./",
            "deseq2_path_prefix": "./"
        ],
        COMPARE_GROUPS.out.artifacts.mix(deseq2_results).collect()
    )

    ch_core_atlas_input_files = core_atlas.concat(
        core_atlas_epithelial_cells, core_atlas_tumor_cells
    ).collect()
    CELL_TYPE_MARKERS_CORE_ATLAS(
        Channel.value(
            [[id: '94a_cell_type_markers_core_atlas'], file("${baseDir}/analyses/90_plots_and_comparisons/94a_cell_type_markers_core_atlas.py")]
        ),
        ch_core_atlas_input_files.map{ core, core_epi, core_tumor -> [
            "main_adata": core.name,
            "epithelial_adata": core_epi.name,
            "tumor_adata": core_tumor.name
        ]},
        ch_core_atlas_input_files
    )
    OVERVIEW_PLOTS_CORE_ATLAS(
        Channel.value(
            [[id: '94b_overview_plots_core_atlas'], file("${baseDir}/analyses/90_plots_and_comparisons/94b_overview_plots_core_atlas.py")]
        ),
        ch_core_atlas_input_files.map{ core, core_epi, core_tumor -> [
            "main_adata": core.name,
            "epithelial_adata": core_epi.name,
            "tumor_adata": core_tumor.name
        ]},
        ch_core_atlas_input_files
    )

    ch_extended_atlas_input_files = extended_atlas.concat(
        // yes, this does need these two files from the *core* atlas. Check notebook for more details.
        core_atlas_epithelial_cells, core_atlas_tumor_cells
    ).collect()
    OVERVIEW_PLOTS_EXTENDED_ATLAS(
        Channel.value(
            [[id: '94c_overview_plots_extended_atlas'], file("${baseDir}/analyses/90_plots_and_comparisons/94c_overview_plots_extended_atlas.py")]
        ),
        ch_extended_atlas_input_files.map{ extended, core_epi, core_tumor -> [
            "main_adata": extended.name,
            "epithelial_adata": core_epi.name,
            "tumor_adata": core_tumor.name
        ]},
        ch_extended_atlas_input_files
    )

    CELL_TYPE_SIGNATURES(
       Channel.value(
            [[id: '94d_cell_type_signatures'], file("${baseDir}/analyses/90_plots_and_comparisons/94d_cell_type_signatures.py")]
        ),
        extended_atlas.map { extended -> [
            "adata_path": extended.name
        ]},
        extended_atlas
    )
    ch_cell_type_signatures = CELL_TYPE_SIGNATURES.out.artifacts.flatten().filter{ it -> it.name.equals("cell_type_major_signatures.csv") }

    COMPARE_PLATFORMS(
        Channel.value(
            [[id: 'compare_platforms'], file("${baseDir}/analyses/90_plots_and_comparisons/95_compare_platforms.py")]
        ),
        extended_atlas.map{ f -> [
            "main_adata": f.name
        ]},
        extended_atlas
    )

    ch_response_to_ici_input_files = ch_neutro_sigs.concat(
        ch_cell_type_signatures,
        Channel.fromPath("${baseDir}/data/14_ici_treatment/Genentech_for_scissor/genentech.rds"),
        Channel.fromPath("${baseDir}/data/14_ici_treatment/Genentech_for_scissor/genentech_clinical_data.tsv")
    ).collect()
    RESPONSE_TO_ICI(
        Channel.value(
            [[id: 'response_to_ici'], file("${baseDir}/analyses/90_plots_and_comparisons/97_response_to_icb.Rmd")]
        ),
        ch_response_to_ici_input_files.map{ sigs_neutro, sigs_cell_type, tpm, meta -> [
            "neutro_sigs": sigs_neutro.name,
            "cell_type_major_sigs": sigs_cell_type.name,
            "ici_tpm": tpm.name,
            "ici_meta": meta.name
        ]},
        ch_response_to_ici_input_files
    )

    ch_cpdb_analysis_input_files = extended_atlas.concat(
        adata_neutrophil_clusters,
        Channel.fromPath("${baseDir}/tables/cellphonedb_2022-04-06.tsv"),
        deseq2_results_neutro_clusters,
        deseq2_results
    ).collect()
    CPDB_ANALYSIS(
        Channel.value(
            [[id: 'cell2cell'], file("${baseDir}/analyses/90_plots_and_comparisons/99_cpdb_analysis.py")]
        ),
        ch_cpdb_analysis_input_files.map{ it -> [
            "main_adata": it[0].name,
            "adata_n": it[1].name,
            "path_cpdb": it[2].name,
            "deseq2_path_neutro_clusters": it[3].name,
            "deseq2_path_prefix": "./"
        ]},
        ch_cpdb_analysis_input_files
    )
    CPDB_CIRCOS_PLOT(
         Channel.value(
            [[id: 'cell2cell_circos'], file("${baseDir}/analyses/90_plots_and_comparisons/99b_cpdb_circosplot.Rmd")]
        ),
        [],
        CPDB_ANALYSIS.out.artifacts
    )


    SCCODA_CONDITION(
        Channel.value([
            [id: "sccoda_condition"],
            file("${baseDir}/analyses/90_plots_and_comparisons/98a_cell_type_composition_luad_lusc.py")
        ]),
        [
            "cell_type_column": "cell_type_major",
            "reference_cell_type": "Tumor cells",
            "mcmc_iterations": 500000,
            "main_adata": "full_atlas_merged.h5ad"
        ],
        extended_atlas.collect()
    )


}

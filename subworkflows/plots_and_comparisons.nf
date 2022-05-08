include {
    JUPYTERNOTEBOOK as NEUTROPHIL_ANALYSIS;
    JUPYTERNOTEBOOK as NEUTROPHIL_ANALYSIS_VELOCYTO;
    JUPYTERNOTEBOOK as STRATIFY_PATIENTS_FIGURES;
    JUPYTERNOTEBOOK as COMPARE_GROUPS;
    JUPYTERNOTEBOOK as COMPARE_GROUPS_PLOTS;
    JUPYTERNOTEBOOK as SCCODA_CONDITION;
    JUPYTERNOTEBOOK as SCCODA_ORIGIN
} from "../modules/local/jupyternotebook/main.nf"

workflow plots_and_comparisons {
    take:
    adata_annotated
    adata_neutrophil_clusters
    patient_stratification

    main:

    ch_neutrophil_analysis_input_files = adata_annotated.concat(
        adata_neutrophil_clusters,
        patient_stratification,
        Channel.fromPath("${baseDir}/tables/gene_annotations/neutro_phenotype_genesets.xlsx"),
        Channel.fromPath("${baseDir}/tables/gene_annotations/neutro_recruitment_chemokines.xlsx")
    ).collect()
    NEUTROPHIL_ANALYSIS(
        Channel.value(
            [[id: 'neutrophil_analysis'], file("${baseDir}/analyses/90_plots_and_comparisons/96a_neutrophils.py")]
        ),
        ch_neutrophil_analysis_input_files.map { adata, adata_n, patient_strat, genesets1, genesets2 -> [
            "adata_n_path": adata_n.name,
            "adata_path": adata.name,
            "patient_stratification_path": patient_strat.name,
            "neutro_geneset_path": genesets1.name,
            "neutro_recruitment_geneset_path": genesets2.name
        ]},
        ch_neutrophil_analysis_input_files
    )

    ch_velocyto_input_files = adata_neutrophil_clusters.concat(
        Channel.fromPath("${baseDir}/data/11_own_datasets/velocyto/")
    ).view().collect()
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
        adata_annotated.mix(patient_stratification).collect()
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
        adata_annotated.collect()
    )


}

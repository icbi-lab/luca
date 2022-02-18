include { JUPYTERNOTEBOOK as COMPARE_GROUPS } from "../modules/local/jupyternotebook/main.nf"

workflow plots_and_comparisons {
    take:
    adata_annotated
    adata_cpdb
    patient_stratification

    main:
    COMPARE_GROUPS(
        [[id: 'compare_groups'], file("${baseDir}/analyses/90_plots_and_comparisons/91_compare_groups.py")],
        Channel.from([
            "tumor_normal",
            "infiltration_status",
            "infiltration_type",
            "patient_immune_infiltration",
            "patient_immune_infiltration_condition",
            "patient_immune_infiltration_treatment_coding",
            "patient_immune_infiltration_treatment_coding_condition",
            "luad_lscc",
            "early_advanced"
        ]).map {
            it -> [
                "comparison": it,
                "adata_in": "full_atlas_merged.h5ad",
                "adata_in_cpdb": "adata_cpdb.h5ad",
                "stratification_csv": "patient_stratification.csv"
            ]
        },
        adata_annotated.mix(adata_cpdb, patient_stratification).collect()
    )

}

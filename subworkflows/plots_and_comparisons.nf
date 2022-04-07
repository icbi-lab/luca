include { JUPYTERNOTEBOOK as COMPARE_GROUPS } from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as SCCODA_CONDITION; JUPYTERNOTEBOOK as SCCODA_ORIGIN } from "../modules/local/jupyternotebook/main.nf"

workflow plots_and_comparisons {
    take:
    adata_annotated
    patient_stratification

    main:
    COMPARE_GROUPS(
        [[id: 'compare_groups'], file("${baseDir}/analyses/90_plots_and_comparisons/91_compare_groups.py")],
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

    SCCODA_ORIGIN(
        [
            [id: "sccoda_condition"],
            file("${baseDir}/analyses/90_plots_and_comparisons/98a_cell_type_composition_luad_lusc.py")
        ],
        [
            "cell_type_column": "cell_type_major",
            "reference_cell_type": "Tumor cells",
            "mcmc_iterations": 500000,
            "main_adata": "full_atlas_merged.h5ad"
        ],
        adata_annotated.collect()
    )
    SCCODA_CONDITION(
        [
            [id: "sccoda_origin"],
            file("${baseDir}/analyses/90_plots_and_comparisons/98b_cell_type_composition_tumor_normal.py")
        ],
        [
            "cell_type_column": "cell_type_major",
            "reference_cell_type": "Stromal",
            "mcmc_iterations": 500000,
            "main_adata": "full_atlas_merged.h5ad"
        ],
        adata_annotated.collect()
    )

}

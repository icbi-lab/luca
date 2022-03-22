include { JUPYTERNOTEBOOK as COMPARE_GROUPS } from "../modules/local/jupyternotebook/main.nf"
include { JUPYTERNOTEBOOK as SCCODA } from "../modules/local/jupyternotebook/main.nf"

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
            "patient_immune_infiltration",
            "patient_immune_infiltration_treatment_coding",
            "patient_immune_infiltration_treatment_coding_condition",
            "patient_immune_infiltration_treatment_coding_condition_random",
            "luad_lusc",
            "early_advanced",
            "early_advanced_condition",
        ]).map {
            it -> [
                "comparison": it,
                "adata_in": "full_atlas_neutrophil_clusters.h5ad",
                "adata_in_cpdb": "adata_cpdb.h5ad",
                "stratification_csv": "patient_stratification.csv"
            ]
        },
        adata_annotated.mix(adata_cpdb, patient_stratification).collect()
    )

    ch_sccoda_params = adata_annotated.map{ it -> it.name }.combine([500000]).combine(
        ["cell_type_major", "cell_type_neutro"]
    ).combine(["Tumor cells", "Stromal"]).map{
        adata, mcmc, col, ref -> [
            "cell_type_column": col,
            "reference_cell_type": ref,
            "mcmc_iterations": mcmc,
            "main_adata": adata
        ]
    }

    SCCODA(
        ch_sccoda_params.map{params -> [
            [id: "${params['cell_type_column']}_${params['reference_cell_type']}"],
            file("${baseDir}/analyses/90_plots_and_comparisons/98_cell_type_composition.py")
        ]},
        ch_sccoda_params,
        adata_annotated.collect()
    )

}

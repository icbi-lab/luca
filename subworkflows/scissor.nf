include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { H5AD_TO_SCE }  from "../modules/local/scconversion/main.nf"
include { SCISSOR as SCISSOR_TCGA; SCISSOR as SCISSOR_GENENTECH } from "../modules/local/scissor.nf"


workflow scissor {
    take: adata_annotated

    main:
    ch_adata_integrated = adata_annotated.map{ it -> [it.baseName, it]}
    SPLIT_ANNDATA(ch_adata_integrated, "patient")
    H5AD_TO_SCE(SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it]})

    ch_sce = H5AD_TO_SCE.out.sce
    // SCISSOR_TCGA(
    //     ch_sce,
    //     file("$baseDir/data/13_tcga/for_scissor/nsclc_primary_tumor.rds", checkIfExists: true),
    //     file("$baseDir/tables/tcga/clinical_data_for_scissor.tsv", checkIfExists: true),
    //     "TCGA_patient_barcode",
    //     Channel.from(
    //         [
    //             "--column tumor_stage",
    //             "--column response_to_chemotherapy",
    //             // "--column kras_mutation",
    //             // "--column braf_mutation",
    //             // "--column egfr_mutation",
    //             "--column tumor_type",
    //             "--surv_time time --surv_status status",

    //             "--tumor_type LUAD --column response_to_chemotherapy",
    //             "--tumor_type LUAD --column braf_mutation",
    //             "--tumor_type LUAD --column kras_mutation", // only enough patients in LUAD
    //             "--tumor_type LUAD --column egfr_mutation",
    //             "--tumor_type LUAD --column tp53_mutation",
    //             "--tumor_type LUAD --column stk11_mutation",
    //             "--tumor_type LUAD --column stk11_kras_mutation", // only enough patients in LUAD
    //             "--tumor_type LUAD --column tumor_stage",
    //             "--tumor_type LUAD --column random",
    //             "--tumor_type LUAD --surv_time time --surv_status status",

    //             "--tumor_type LUSC --column response_to_chemotherapy",
    //             "--tumor_type LUSC --column braf_mutation",
    //             "--tumor_type LUSC --column egfr_mutation",
    //             "--tumor_type LUSC --column tp53_mutation",
    //             "--tumor_type LUSC --column stk11_mutation",
    //             "--tumor_type LUSC --column tumor_stage",
    //             "--tumor_type LUSC --column random",
    //             "--tumor_type LUSC --surv_time time --surv_status status"
    //         ]
    //     )
    // )

    SCISSOR_GENENTECH(
        ch_sce,
        file("$baseDir/data/14_ici_treatment/Genentech_for_scissor/genentech.rds", checkIfExists: true),
        file("$baseDir/data/14_ici_treatment/Genentech_for_scissor/genentech_clinical_data.tsv", checkIfExists: true),
        "sample_id",
        Channel.from(
            [
                "--column response_to_chemotherapy",
                "--column response_to_ici",
                "--surv_time time --surv_status status",
                "--surv_time time_ici --surv_status status_ici",
                "--surv_time time_chemo --surv_status status_chemo",

                "--tumor_type LUAD --column response_to_ici",
                "--tumor_type LUAD --column response_to_chemotherapy",
                "--tumor_type LUAD --surv_time time --surv_status status",
                "--tumor_type LUAD --surv_time time_ici --surv_status status_ici",
                "--tumor_type LUAD --surv_time time_chemo --surv_status status_chemo",

                "--tumor_type LUSC --column response_to_ici",
                "--tumor_type LUSC --column response_to_chemotherapy",
                "--tumor_type LUSC --surv_time time --surv_status status",
                "--tumor_type LUSC --surv_time time_ici --surv_status status_ici",
                "--tumor_type LUSC --surv_time time_chemo --surv_status status_chemo"
            ]
        )
    )
}

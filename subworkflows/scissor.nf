include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { FILTER_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { H5AD_TO_SCE }  from "../modules/local/scconversion/main.nf"
include { SCISSOR as SCISSOR_TCGA } from "../modules/local/scissor.nf"
include {
    RMARKDOWNNOTEBOOK as VALIDATE_DECONVOLUTION;
} from "../modules/local/rmarkdownnotebook/main.nf"
include {
    JUPYTERNOTEBOOK as SCISSOR_ANALYSIS;
} from "../modules/local/jupyternotebook/main.nf"

workflow scissor {
    take: adata_annotated

    main:
    ch_adata_integrated = adata_annotated.map{ it -> [it.baseName, it]}
    FILTER_ANNDATA(
        ch_adata_integrated,
        """lambda x: (x['origin'] == 'tumor_primary')"""
    )
    SPLIT_ANNDATA(FILTER_ANNDATA.out.adata, "patient")
    H5AD_TO_SCE(SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it]})

    ch_sce = H5AD_TO_SCE.out.sce
    SCISSOR_TCGA(
        ch_sce,
        file("$baseDir/data/13_tcga/for_scissor/nsclc_primary_tumor.rds", checkIfExists: true),
        file("$baseDir/tables/tcga/clinical_data_for_scissor.tsv", checkIfExists: true),
        "TCGA_patient_barcode",
        Channel.from(
            [
                // "--column kras_mutation",
                // "--column braf_mutation",
                // "--column egfr_mutation",

                "--column tumor_type",
                "--surv_time time --surv_status status",

                "--tumor_type LUAD --column braf_mutation",
                "--tumor_type LUAD --column kras_mutation", // only enough patients in LUAD
                "--tumor_type LUAD --column egfr_mutation",
                "--tumor_type LUAD --column tp53_mutation",
                "--tumor_type LUAD --column stk11_mutation",
                "--tumor_type LUAD --column random",
                "--tumor_type LUAD --surv_time time --surv_status status",

                "--tumor_type LUSC --column braf_mutation",
                "--tumor_type LUSC --column egfr_mutation",
                "--tumor_type LUSC --column tp53_mutation",
                "--tumor_type LUSC --column stk11_mutation",
                "--tumor_type LUSC --column random",
                "--tumor_type LUSC --surv_time time --surv_status status"
            ]
        )
    )

    ch_scissor_analysis_input_files = adata_annotated.concat(
        Channel.fromPath("$baseDir/tables/tcga/clinical_data_for_scissor.tsv"),
        SCISSOR_TCGA.out.scissor_cells
    ).collect()
    SCISSOR_ANALYSIS(
        Channel.value(
            [[id: 'scissor_analysis'], file("${baseDir}/analyses/50_scissor/51_scissor_analysis.py")]
        ),
        ch_scissor_analysis_input_files.map{ it -> [
            "adata_in": it[0].name,
            "path_clinical_data": it[1].name,
            "path_scissor": "./"
        ]},
        ch_scissor_analysis_input_files
    )

    ch_validate_deconvolution_input_files = Channel.fromPath("${baseDir}/tables/tcga/clinical_data_for_scissor.tsv").concat(
        Channel.fromPath("${baseDir}/data/13_tcga/for_scissor/nsclc_primary_tumor.rds"),
        Channel.fromPath("${baseDir}/tables/tcga/mmc1.xlsx"),
    ).collect()
    VALIDATE_DECONVOLUTION(
        Channel.value(
            [[id: 'scissor_validate_deconvolution'], file("${baseDir}/analyses/50_scissor/52_validate_deconvolution.Rmd")]
        ),
        ch_validate_deconvolution_input_files.map{ clinical, tpm, mmc1 -> [
            "clinical_data": clinical.name,
            "tcga_tpm": tpm.name,
            "tcga_meta": mmc1.name
        ]},
        ch_validate_deconvolution_input_files
    )

}

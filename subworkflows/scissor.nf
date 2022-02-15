include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { H5AD_TO_SCE }  from "../modules/local/scconversion/main.nf"

process SCISSOR {
    input:
    tuple val(id), path(sce)
    path(bulk_tpm)
    path(metadata)
    each options

    output:
    path "scissor*.tsv", emit: scissor_cells, optional: true
    path "*.log", emit: log

    script:
    // using this instead of `errorStrategy` in order to also cache failed processes
    // (they will always fail due to characteristics of the data, e.g. too few cells)
    ignore_exit_code = task.ext.ignore_error ? "|| true" : ""
    """
    scissor_single_sample.R --bulk_tpm $bulk_tpm --sce $sce --metadata $metadata \\
        --sample_col=TCGA_patient_barcode \\
        $options > ${id}_${options.replace("-","").replace(" ", "_")}.log 2>&1 $ignore_exit_code
    """
}

workflow scissor {
    take: adata_annotated

    main:
    ch_adata_integrated = adata_annotated.map{ it -> [it.baseName, it]}
    SPLIT_ANNDATA(ch_adata_integrated, "patient")
    H5AD_TO_SCE(SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it]})

    ch_sce = H5AD_TO_SCE.out.sce
    SCISSOR(
        ch_sce,
        file("$baseDir/data/13_tcga/for_scissor/nsclc_primary_tumor.rds", checkIfExists: true),
        file("$baseDir/tables/tcga/clinical_data_for_scissor.tsv", checkIfExists: true),
        Channel.from(
            [
                // "--column tumor_stage",
                // "--column kras_mutation",
                // "--column braf_mutation",
                // "--column egfr_mutation",
                "--column tumor_type",
                // "--surv_time time --surv_status status",
                "--tumor_type LUAD --column kras_mutation",
                "--tumor_type LUAD --column braf_mutation",
                "--tumor_type LUAD --column egfr_mutation",
                "--tumor_type LUAD --column tumor_stage",
                "--tumor_type LUAD --column tp53_mutation",
                "--tumor_type LUAD --column random",
                "--tumor_type LUAD --surv_time time --surv_status status",
                "--tumor_type LUSC --column kras_mutation",
                "--tumor_type LUSC --column braf_mutation",
                "--tumor_type LUSC --column egfr_mutation",
                "--tumor_type LUSC --column tp53_mutation",
                "--tumor_type LUSC --column tumor_stage",
                "--tumor_type LUSC --column random",
                "--tumor_type LUSC --surv_time time --surv_status status"
            ]
        )
    )
}

include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { H5AD_TO_SCE }  from "../modules/local/scconversion/main.nf"
include { JUPYTERNOTEBOOK as ITH_ANALYSIS } from "../modules/local/jupyternotebook/main.nf"

process RUN_SCEVAN {
    input:
    tuple val(id), path(input_file)

    output:
    path "${id}", emit: res_dir
    path "${id}.log", emit: log

    script:
    // using this instead of `errorStrategy` in order to also cache failed processes
    // (they will always fail due to characteristics of the data, e.g. too few cells)
    ignore_exit_code = task.ext.ignore_error ? "|| true" : ""
    """
    scevan_parallel.R ${input_file} ${task.cpus} ${input_file.baseName} > ${id}.log 2>&1 $ignore_exit_code
    mv scevan_result.csv output
    mv output ${id}
    """
}

process RUN_INFERCNVPY {
    input:
    tuple val(id), path(input_file)
    path(gtffile)

    output:
    path "*", emit: output_file

    script:
    // using this instead of `errorStrategy` in order to also cache failed processes
    // (they will always fail due to characteristics of the data, e.g. too few cells)
    ignore_exit_code = task.ext.ignore_error ? "|| true" : ""
    """
    infercnvpy_parallel.py ${input_file} ${gtffile} > ${id}.log 2>&1 $ignore_exit_code
    """
}

workflow infercnv {
    take:
    adata_annotated
    patient_stratification_table

    main:
    ch_adata_annotated = adata_annotated.map{ it -> [it.baseName, it]}
    SPLIT_ANNDATA(ch_adata_annotated, "patient")

    ch_adatas_by_patient = SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it]}
    H5AD_TO_SCE(ch_adatas_by_patient)

    // RUN_INFERCNVPY(
    //     ch_adatas_by_patient,
    //     file("${baseDir}/data/10_references/gencode.v33.primary_assembly.annotation.gtf", checkIfExists: true)
    // )
    RUN_SCEVAN(H5AD_TO_SCE.out.sce)

    ch_ith_analysis_input_files = adata_annotated.concat(
        patient_stratification_table,
        RUN_SCEVAN.out.res_dir
    ).collect()
    ITH_ANALYSIS(
        Channel.value(
            [[id: 'ith_analysis'], file("${baseDir}/analyses/60_cnv_analysis/ith_analysis.py")]
        ),
        ch_ith_analysis_input_files.map{ it -> [
            "adata_in": it[0].name,
            "stratification_csv": it[1].name,
            "path_scevan": "./"
        ]},
        ch_ith_analysis_input_files
    )
}



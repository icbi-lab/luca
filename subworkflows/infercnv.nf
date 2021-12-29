include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { H5AD_TO_SCE }  from "../modules/local/scconversion/main.nf"

process RUN_SCEVAN {
    input:
    tuple val(id), path(input_file)

    output:
    path "*", emit: output_file

    errorStrategy 'ignore'

    script:
    """
    scevan_parallel.R ${input_file} ${task.cpus} ${input_file.baseName}
    """
}

// process run_copykat {
//     input:
//     path input_file

//     output:
//     path "*", emit: output_file

//     errorStrategy 'ignore'

//     script:
//     """
//     copykat_parallel.R ${input_file} ${task.cpus} ${input_file.baseName}
//     """
// }

process RUN_INFERCNVPY {
    input:
    tuple val(id), path(input_file)

    output:
    path "*", emit: output_file

    errorStrategy 'ignore'

    script:
    """
    infercnvpy_parallel.py ${input_file}
    """
}

workflow infercnv {
    take: adata_annotated

    main:
    ch_adata_annotated = Channel.value([adata_annotated.baseName, adata_annotated])
    SPLIT_ANNDATA(ch_adata_annotated, "patient")

    ch_adatas_by_patient = SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it]}
    H5AD_TO_SCE(ch_adatas_by_patient)

    RUN_INFERCNVPY(ch_adatas_by_patient)
    RUN_SCEVAN(H5AD_TO_SCE.out.sce)
}



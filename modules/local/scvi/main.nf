#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process SCVI {
    input:
        tuple val(id), path(input_adata)
        each use_highly_variable
        tuple val(batch_key), val(hvg_batch_key), val(labels_key)

    output:
        tuple val(id), path("*_integrated_scvi.h5ad"), emit: adata
        tuple val(id), path("*_scvi_model"), emit: scvi_model

    script:
    def suffix = use_highly_variable == 0 ? "all_genes" : "hvg"
    def labels_key_arg = labels_key ? "--labels_key ${labels_key}" : ""
    """
    export CUDA_VISIBLE_DEVICES=\$((0 + \$RANDOM % 2))
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}  \\
        MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus}  \\
        MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus}
    integrate_scvi.py \\
        ${input_adata} \\
        --adata_out ${id}_${suffix}_integrated_scvi.h5ad \\
        --model_out ${id}_${suffix}_scvi_model \\
        --use_hvg ${use_highly_variable} \\
        --hvg_batch_key ${hvg_batch_key} \\
        --batch_key ${batch_key} \\
        ${labels_key_arg}
    """
}

process SCANVI {
    input:
        tuple val(id), path(input_adata), path(input_model)
        val batch_key
        val labels_key

    output:
        tuple val(id), path("*_integrated_scanvi.h5ad"), emit: adata
        tuple val(id), path("*scanvi_model"), emit: scvi_model

    script:
    """
    export CUDA_VISIBLE_DEVICES=\$((0 + \$RANDOM % 2))
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}  \\
        MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus}  \\
        MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus}
    integrate_scanvi.py \\
        ${input_adata} \\
        ${input_model} \\
        --adata_out ${input_adata.baseName}_integrated_scanvi.h5ad \\
        --model_out ${input_adata.baseName}_scanvi_model \\
        --batch_key ${batch_key} \\
        --labels_key ${labels_key}
    """
}


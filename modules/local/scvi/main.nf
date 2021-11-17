#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { initOptions; saveFiles; getSoftwareName } from './functions'

process SCVI {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cpus 4
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    // container "containers/sc-integrate2.sif"
    // // support for nvidia https://lucacozzuto.medium.com/using-gpu-within-nextflow-19cd185d5e69
    // containerOptions = "--nv"
    clusterOptions '-V -S /bin/bash -q all.q@apollo-15'

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
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // don't run on GPU due to pytorch internal runtime error
    // runs in reasonable time (20min-ish) without a GPU on a
    // full compute node
    cpus 44
    // container "containers/sc-integrate2.sif"
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"

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


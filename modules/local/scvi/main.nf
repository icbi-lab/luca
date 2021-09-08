#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { initOptions; saveFiles; getSoftwareName } from './functions'

process SCVI {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cpus 4
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    clusterOptions '-V -S /bin/bash -q all.q@apollo-15'

    input:
        path input_adata
        each use_highly_variable
        val batch_key

    output:
        path("integrated*.h5ad"), emit: adata
        path("scvi_model*"), emit: scvi_model

    script:
    def suffix = use_highly_variable == 0 ? "all_genes" : "hvg"
    """
    export CUDA_VISIBLE_DEVICES=\$((0 + \$RANDOM % 2))
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}  \\
        MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus}  \\
        MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus}
    integrate_scanvi.py \\
        ${input_adata} \\
        --adata_out integrated_${input_adata.baseName}_${suffix}.h5ad \\
        --model_out scvi_model_${input_adata.baseName}_${suffix} \\
        --use_hvg ${use_highly_variable} \\
        --batch_key ${batch_key}
    """
}


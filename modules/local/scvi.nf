#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SCVI {
    publishDir "../../data/50_integrate_scrnaseq_data/52_run_scanvi", mode: "copy"

    cpus 4
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    clusterOptions '-V -S /bin/bash -q all.q@apollo-15'

    input:
        path input_adata
        each use_highly_variable

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
        integrated_${input_adata.baseName}_${suffix}.h5ad \\
        scvi_model_${input_adata.baseName}_${suffix} \\
        ${use_highly_variable}
    """
}


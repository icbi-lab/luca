#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SCANVI {
    publishDir "../../data/50_integrate_scrnaseq_data/"

    cpus 4
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    clusterOptions '-V -S /bin/bash -q all.q@apollo-15'

    input:
        path input_adata

    output:
        path("*.integrated.h5ad")

    script:
    """
    export CUDA_VISIBLE_DEVICES=\$((0 + \$RANDOM % 2))
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}  \\
        MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus}  \\
        MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus}
    integrate_scanvi.py ${input_adata} ${input_adata.baseName}.integrated.h5ad 
    """
}

workflow {
    SCANVI(
        file("../../data/50_integrate_scrnaseq_data/merged*.h5ad"),
    )
}

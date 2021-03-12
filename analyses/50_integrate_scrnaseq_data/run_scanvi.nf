#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SCANVI {
    publishDir "../../data/50_integrate_scrnaseq_data/", mode: "copy"

    cpus 4
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    clusterOptions '-V -S /bin/bash -q all.q@apollo-15'

    input:
        path input_adata

    output:
        path("integrated*.h5ad")
        path("scvi_model*")

    script:
    """
    export CUDA_VISIBLE_DEVICES=\$((0 + \$RANDOM % 2))
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}  \\
        MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus}  \\
        MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus}
    integrate_scanvi.py ${input_adata} integrated_${input_adata.baseName}.h5ad scvi_model_${input_adata.baseName}
    """
}

workflow {
    SCANVI(
        Channel.fromPath("../../data/50_integrate_scrnaseq_data/10_merge_all/merged*.h5ad"),
    )
}

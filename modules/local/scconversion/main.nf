#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { initOptions; saveFiles; getSoftwareName } from './functions'

process H5AD_TO_SEURAT {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cpus 1
    conda "/data/scratch/sturm/conda/envs/2020-pircher-seuratdisk"
    stageInMode 'link'

    input:
        tuple val(id), path(input_adata)

    output:
        tuple val(id), path("*.h5seurat"), emit: h5seurat

    script:
    """
    #!/usr/bin/env Rscript
    SeuratDisk::Convert("${input_adata}", dest="h5seurat")
    """
}


process SEURAT_TO_SCE {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cpus 1
    conda "/data/scratch/sturm/conda/envs/2020-pircher-seuratdisk"

    input:
        tuple val(id), path(input_seurat)

    output:
        tuple val(id), path("*.rds"), emit: sce

    // TODO: the conversion fails for some objects, likely this is related to `.uns`
    script:
    """
    #!/usr/bin/env Rscript

    library(Seurat)
    seurat_obj = SeuratDisk::LoadH5Seurat("${input_seurat}")
    sce_obj = as.SingleCellExperiment(seurat_obj)
    saveRDS(sce_obj, file="${id}.rds", compress=FALSE)
    """
}

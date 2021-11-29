#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process H5AD_TO_SEURAT {
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

process H5AD_TO_SCE {
    // requires rpy2 <= 3.4.2, see https://github.com/theislab/anndata2ri/issues/63
    conda "/data/scratch/sturm/conda/envs/2020-pircher-seuratdisk"
    stageInMode 'link'

    input:
        tuple val(id), path(input_adata)

    output:
        tuple val(id), path("*.rds"), emit: sce

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    adata = sc.read_h5ad("${input_adata}")

    try:
        from rpy2.robjects.packages import importr
        from rpy2.robjects.conversion import localconverter
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2 import robjects as ro
        import anndata2ri
        r_base = importr("base")
    except ImportError:
        raise ImportError("rpy2 and anndata2ri need to be installed. ")

    with localconverter(anndata2ri.converter):
        sce = ro.conversion.py2rpy(adata)

    r_base.saveRDS(sce, file="${id}.rds", compress=False)
    """
}


process SEURAT_TO_SCE {
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
    # misc == adata.uns; Excluding it because it triggered an error when converting to SCE
    seurat_obj = SeuratDisk::LoadH5Seurat("${input_seurat}", misc=FALSE)
    sce_obj = as.SingleCellExperiment(seurat_obj)
    saveRDS(sce_obj, file="${id}.rds", compress=FALSE)
    """
}


process SPLIT_ANNDATA {
    input:
        tuple val(id), path(input_adata)
        val(split_variable)

    output:
        path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    import re

    adata = sc.read_h5ad("${input_adata}")
    for v in adata.obs["${split_variable}"].unique():
        tmp_adata = adata[adata.obs["${split_variable}"] == v, :].copy()
        v_sanitized = re.sub("[^a-z0-9_-]", "_", v.lower())
        tmp_adata.write_h5ad(f"./${id}_{v_sanitized}.h5ad")
    """
}

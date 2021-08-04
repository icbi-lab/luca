
include { initOptions; saveFiles; getSoftwareName } from './functions'


process NEIGHBORS {
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    cpus 8

    input:
    path adata
    val use_rep

    output:
    path "*.h5ad", emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.pp.neighbors(adata, use_rep="${use_rep}")
    adata.write_h5ad("${adata.baseName}.neighbors.h5ad")
    """
}

process UMAP {
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    cpus 8

    input:
    path adata

    output:
    path "*.h5ad", emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.umap(adata)
    adata.write_h5ad("${adata.simpleName}.umap.h5ad")
    """
}

process LEIDEN {
    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    cpus 1

    input:
    path adata
    val resolution

    output:
    path "*.h5ad", emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.leiden(adata, resolution=${resolution})
    adata.write_h5ad("${adata.simpleName}.leiden.h5ad")
    """
}

process MERGE_UMAP_LEIDEN {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    cpus 1

    input:
    path adata_umap
    path adata_leiden

    output:
    path "*.h5ad", emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata_umap = sc.read_h5ad("${adata_umap}")
    adata_leiden = sc.read_h5ad("${adata_leiden}")
    adata_umap.obs["leiden"] = adata_leiden.obs["leiden"]
    adata_umap.write_h5ad("${adata_umap.simpleName}.umap_leiden.h5ad")
    """
}



workflow NEIGHBORS_LEIDEN_UMAP {
    take:
    adata
    neihbors_rep
    leiden_res

    main:
    NEIGHBORS(adata, neihbors_rep)
    UMAP(NEIGHBORS.out.adata)
    LEIDEN(NEIGHBORS.out.adata, leiden_res)
    MERGE_UMAP_LEIDEN(UMAP.out.adata, LEIDEN.out.adata)

    emit:
    adata = MERGE_UMAP_LEIDEN.out.adata
}

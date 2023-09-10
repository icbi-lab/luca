

process NEIGHBORS {
    container = "biomedbigdata/sc-integration"
    cpus 8

    input:
    tuple val(id), path(adata)
    val use_rep

    output:
    tuple val(id), path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.pp.neighbors(adata, use_rep="${use_rep}")
    adata.write_h5ad("${id}.neighbors.h5ad")
    """
}

process UMAP {
    container = "biomedbigdata/sc-integration"
    cpus 8

    input:
    tuple val(id), path(adata)

    output:
    tuple val(id), path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.umap(adata)
    adata.write_h5ad("${id}.umap.h5ad")
    """
}

process LEIDEN {
    container = "biomedbigdata/sc-integration"
    cpus 1

    input:
    tuple val(id), path(adata)
    each resolution

    output:
    tuple val(id), val(resolution), path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    adata = sc.read_h5ad("${adata}")
    sc.tl.leiden(adata, resolution=${resolution})
    adata.write_h5ad("${id}.res_${resolution}.leiden.h5ad")
    """
}

process MERGE_UMAP_LEIDEN {
    container = "biomedbigdata/sc-integration"
    cpus 1

    input:
    tuple val(id), path(adata_umap), val(leiden_resolutions), path(adata_leiden)

    output:
    tuple val(id), path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    from threadpoolctl import threadpool_limits

    threadpool_limits(${task.cpus})
    sc.settings.n_jobs = ${task.cpus}

    resolutions = ${leiden_resolutions}
    if not isinstance(resolutions, list):
        resolutions = [resolutions]
    leiden_adatas = "${adata_leiden}".split(" ")

    adata_umap = sc.read_h5ad("${adata_umap}")
    for res, adata_path in zip(resolutions, leiden_adatas):
        tmp_adata = sc.read_h5ad(adata_path)
        adata_umap.obs[f"leiden_{res:.2f}"] = tmp_adata.obs["leiden"]
    adata_umap.write_h5ad("${id}.umap_leiden.h5ad")
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

    MERGE_UMAP_LEIDEN(UMAP.out.adata.join(LEIDEN.out.groupTuple()))

    emit:
    adata = MERGE_UMAP_LEIDEN.out.adata
}

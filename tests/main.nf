#!/usr/bin/env nextflow

nextflow.enable.dsl= 2

params.publish_dir_mode = "symlink"
params.outdir = "outdir"
params.options = [:]

include { NEIGHBORS_LEIDEN_UMAP } from "../subworkflows/neighbors_leiden_umap/main.nf"

workflow test_neighbors_leiden_umap_single_res {
    NEIGHBORS_LEIDEN_UMAP(Channel.fromPath('./testdata/adata_all_celltypes.h5ad'), "X_scVI", 0.5)
}

workflow test_neighbors_leiden_umap_multi_res {
    NEIGHBORS_LEIDEN_UMAP(Channel.fromPath('./testdata/adata_all_celltypes.h5ad'), "X_scVI", Channel.from(0.5, 1.0))
}

include { nxfvars } from "./nxfvars.nf"

process SCQC {
    tag = { meta.id }

    input:
    tuple path(scqc_notebook), path(scqc_lib)
    tuple val(meta), path(input_adata)

    output:
    tuple val(meta.id), path(output_adata), emit: adata
    path(output_stats), emit: qc_stats
    path("*.html"), emit: notebook

    script:
    output_adata = "${meta.id}.qc.h5ad"
    output_stats = "${meta.id}.qc_stats.tsv"
    dataset_id = meta.id
    min_counts = meta.min_counts
    max_counts = meta.max_counts
    min_genes = meta.min_genes
    max_genes = meta.max_genes
    max_pct_mito = meta.max_pct_mito
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache

    ${nxfvars(task)}

    nxfvars execute ${scqc_notebook} ${dataset_id}_qc_report.html
    """
}

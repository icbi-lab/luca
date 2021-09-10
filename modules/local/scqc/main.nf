include { nxfvars } from "../nxfvars.nf"
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process SCQC {
    tag = { meta.id }

    publishDir {params.outdir},
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda "/home/sturm/.conda/envs/single-cell-analysis-nf"

    input:
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
    ${nxfvars(task)}

    export PYTHONPATH="${moduleDir}"
    nxfvars execute ${moduleDir}/scqc-notebook.py ${dataset_id}_qc_report.html
    """
}

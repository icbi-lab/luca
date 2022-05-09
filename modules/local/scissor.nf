process SCISSOR {
    input:
    tuple val(id), path(sce)
    path(bulk_tpm)
    path(metadata)
    val(sample_col)
    each options

    output:
    path "scissor*.tsv", emit: scissor_cells, optional: true
    path "*.log", emit: log

    script:
    // using this instead of `errorStrategy` in order to also cache failed processes
    // (they will always fail due to characteristics of the data, e.g. too few cells)
    ignore_exit_code = task.ext.ignore_error ? "|| true" : ""
    """
    scissor_single_sample.R --bulk_tpm $bulk_tpm --sce $sce --metadata $metadata \\
        --sample_col=$sample_col --prefix=$id \\
        $options > ${id}_${options.replace("-","").replace(" ", "_")}.log 2>&1 $ignore_exit_code
    """
}

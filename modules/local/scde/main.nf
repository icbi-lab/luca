#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { initOptions; saveFiles; getSoftwareName } from './functions'


process PREPARE_ANNDATA {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cpus 1
    conda "/data/scratch/sturm/conda/envs/2020-pircher-seuratdisk"

    input:
    tuple val(id), path(input_adata)
    // layer in anndata which contains raw counts. Set to "X" to use `X`.
    val count_layer
    // column containing the dependent group
    val status_col
    // specify which groups to test against which. May be
    //   * a list of strings refering to individual values in `status_col`
    //   * the string "rest" for `control_groups`, referring to all columns
    //     that are not in `case_group`.
    tuple val(case_groups), val(control_groups)

    output:
    // the output anndata contains only X (filled with raw counts), var and obs.
    tuple val(id), path("*.h5ad")


    script:
    def case_groups_py = "[" + case_groups.collect{ "\"$it\"" }.join(", ") + "]"
    def control_groups_py = (control_groups == "rest") ? '"rest"' : "[" + control_groups.collect{ "\"$it\"" }.join(", ") + "]"
    """
    #!/usr/bin/env python

    import scanpy as sc

    layer = "${count_layer}"
    case_groups = ${case_groups_py}
    control_groups = ${control_groups_py}
    status_col = "${status_col}"
    adata = sc.read_h5ad("${input_adata}")

    status_groups = (
        set(list(adata.obs[status_col].unique()))
        if control_groups == "rest"
        else set(case_groups) | set(control_groups)
    )

    # subset to only the categories of interest
    adata = adata[adata.obs[status_col].isin(status_groups), :]

    # create new, reduced anndata with only the information required for DE analysis
    tmp_adata = sc.AnnData(
        X = adata.layers[layer] if layer != "X" else adata.X,
        obs = adata.obs,
        var = adata.var
    )

    tmp_adata.obs[status_col] = [
        "case" if v in case_groups else "control" for v in tmp_adata.obs[status_col]
    ]
    # make categorical
    tmp_adata._sanitize()

    tmp_adata.write_h5ad("${id}_for_de.h5ad")
    """
}


process MAKE_PSEUDOBULK {

    input:
    tuple val(id), path(input_adata)
    // the column containing the biological replicate variable
    val(replicate_column)

    output:
    path("${id}_counts.csv"), emit: counts
    path("${id}_samplesheet.csv"), emit: samplesheet

    script:
    """
    # aggregate counts by sample

    # Create samplesheet with all variables in obs that are unique per group.
    """

}


process DE_EDGER {
    /**
     * For a standard pseudobulk analysis
     */
    input:
    tuple val(id), path(counts), path(samplesheet)
    val(condition_col)
    val(covariate_formula)


    output:
    path("de_res.csv"), emit: de_res

    script:
    """
    # Run edgeR QLF test or LRT Test as descirbed in their vignette
    """
}


process DE_DREAM {
    /**
     * For a pseodubulk analysis with mixed effects
     * See https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/dream.html

     * This is useful if there is a batch effect that is not controlled for
     * by the experimental design
     */

    input:
    tuple val(id), path(counts), path(samplesheet)
    val(condition_col)
    val(covariate_formula)

    output:
    path("de_res.csv"), emit: de_res

    script:
    """
    # Run limma/dream as descirbed in their vignette
    """

}


process DE_MAST_MIXED_EFFECTS {
    /**
     * Mixed effects model with MAST as described by Zimmermann et al.
     */

    input:
    tuple val(id), path(sce)
    val(condition_col)
    val(formula)


    script:
    """
    #!/usr/bin/env Rscript

    options(mc.cores=${task.cpus})
    RhpcBLASctl::blas_set_num_threads(1)
    library(MAST)
    library(readr)

    sce = readRDS(${sce})
    sca = MAST::SceToSingleCellAssay(sce, check_sanity = TRUE)

    res_zlm = zlm(
        # ~ leiden + n_genes_by_counts + (1 | organoid),
        ${formula},
        sca,
        method="glmer",
        ebayes=FALSE,
        strictConvergence=FALSE,
        parallel=TRUE
    )

    contrast = ${condition_col}control
    zlm_summary = MAST::summary(res_zlm, doLRT="")

    summary_dt <- zlm_summary\$datatable
    de_res <- merge(summary_dt[summary_dt\$contrast==contrast
                                  & summary_dt\$component=='logFC', c(1,7,5,6,8)],
                        summary_dt[summary_dt\$contrast==contrast
                                  & summary_dt\$component=='H', c(1,4)],
                        by = 'primerid')

    write_tsv(de_res, "${id}_de_res.tsv")
    """
}


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
    //   * the string "all" for 'case_group', referring to all-against-all
    //     comparisons ("ignores control_groups")
    tuple val(case_groups), val(control_groups)

    output:
    // the output anndata contains only X (filled with raw counts), var and obs.
    tuple val(id), path("*.h5ad"), emit: adata


    script:
    def case_groups_py = (case_groups == "all") ? '"all"' : "[" + case_groups.collect{ "\"$it\"" }.join(", ") + "]"
    def control_groups_py = (control_groups == "rest") ? '"rest"' : "[" + control_groups.collect{ "\"$it\"" }.join(", ") + "]"
    """
    #!/usr/bin/env python

    import scanpy as sc

    layer = "${count_layer}"
    case_groups = ${case_groups_py}
    control_groups = ${control_groups_py}
    status_col = "${status_col}"
    adata = sc.read_h5ad("${input_adata}")

    if case_groups == "all":
        case_groups_iter = [[x] for x in adata.obs[status_col].unique()]
        control_groups = "rest"
    else:
        case_groups_iter = [case_groups]

    for case_groups in case_groups_iter:
        status_groups = (
            set(list(adata.obs[status_col].unique()))
            if control_groups == "rest"
            else set(case_groups) | set(control_groups)
        )

        # subset to only the categories of interest
        adata_sub = adata[adata.obs[status_col].isin(status_groups), :]

        # create new, reduced anndata with only the information required for DE analysis
        tmp_adata = sc.AnnData(
            X = adata_sub.layers[layer] if layer != "X" else adata_sub.X,
            obs = adata_sub.obs,
            var = adata_sub.var
        )

        tmp_adata.obs[status_col] = [
            "case" if v in case_groups else "control" for v in tmp_adata.obs[status_col]
        ]
        # make categorical
        tmp_adata._sanitize()

        tmp_adata.write_h5ad(f"${id}_{'-'.join(case_groups)}_for_de.h5ad")
    """
}


process MAKE_PSEUDOBULK {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cpus 1
    conda "/data/scratch/sturm/conda/envs/2020-pircher-seuratdisk"

    input:
    tuple val(id), path(input_adata)
    // the column containing the biological replicate variable
    val(replicate_col)
    val(condition_col)
    val(min_cells_per_sample)

    output:
    tuple val(id), path("${id}_counts.csv"), path("${id}_samplesheet.csv"), emit: pseudobulk

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc
    import pandas as pd
    import numpy as np

    adata = sc.read_h5ad("${input_adata}")
    replicate_col = "${replicate_col}"
    condition_col = "${condition_col}"
    min_cells_per_sample = $min_cells_per_sample

    # aggregate counts by sample
    bulk_samples = {}
    n_cells = {}
    rep_cond = adata.obs.loc[:, [replicate_col, condition_col]].drop_duplicates()
    for _, replicate, condition in rep_cond.itertuples():
        mask = (adata.obs[replicate_col] == replicate) & (
            adata.obs[condition_col] == condition
        )
        sample_id = f"{replicate}_{condition}"
        n_cells[sample_id] = np.sum(mask)
        bulk_samples[sample_id] = pd.Series(
            np.sum(
                adata[
                    mask,
                    :,
                ].X,
                axis=0,
            ).A1,
            index=adata.var_names,
        )

    bulk_df = pd.DataFrame(bulk_samples)
    bulk_df.index.name = "gene_id"

    # Create samplesheet with all variables in obs that are unique per group.
    keep_cols = np.all(
        adata.obs.groupby([replicate_col, condition_col], observed=True).apply(lambda x: x.nunique()) == 1,
        axis=0
    )
    samplesheet = adata.obs.loc[:, keep_cols].drop_duplicates()
    samplesheet.index = [
        f"{replicate}_{condition}"
        for replicate, condition in
        zip(samplesheet[replicate_col], samplesheet[condition_col])
    ]
    # drop a pre-existing sample column should it exist
    samplesheet.drop(["sample"], axis="columns", inplace=True, errors="ignore")
    samplesheet.index.name = "sample"
    samplesheet["n_cells"] = pd.Series(n_cells)
    samplesheet = samplesheet.loc[samplesheet["n_cells"] > min_cells_per_sample, :]

    # remove all entries that don't have a paired case/control sample
    keep_replicates = samplesheet.groupby([replicate_col]).size() == 2
    keep_replicates = keep_replicates.index[keep_replicates].values
    samplesheet = samplesheet.loc[samplesheet[replicate_col].isin(keep_replicates), :]
    samplesheet.sort_index(inplace=True)
    bulk_df = bulk_df.loc[:, samplesheet.index]

    # Export as CSV
    bulk_df.to_csv("${id}_counts.csv")
    samplesheet.to_csv("${id}_samplesheet.csv")
    """

}


process DE_EDGER {
    /**
     * For a standard pseudobulk analysis
     */
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cpus 1
    conda "/data/scratch/sturm/conda/envs/2020-pircher-seuratdisk"

    input:
    tuple val(id), path(counts), path(samplesheet)
    val(condition_col)
    val(covariate_formula)

    output:
    path("*.tsv"), emit: de_res

    script:
    """
    #!/usr/bin/env Rscript

    options(mc.cores=${task.cpus})
    RhpcBLASctl::blas_set_num_threads(1)

    library(readr)
    library(edgeR)
    library(tibble)

    counts = read_csv("${counts}")
    samplesheet = read_csv("${samplesheet}")

    # Run edgeR QLF test or LRT Test as descirbed in their vignette
    design = model.matrix(~ ${condition_col} ${covariate_formula}, data=samplesheet)
    dge = DGEList(
        counts = column_to_rownames(counts, "gene_id"),
        samples = column_to_rownames(samplesheet, "sample")
    )

    dge = calcNormFactors(dge, design = design)
    dge = estimateDisp(dge, design = design)
    fit = glmQLFit(dge, design = design)

    qlf = glmQLFTest(fit)
    de_res = as_tibble(
        topTags(qlf, n=Inf, adjust.method="BH")\$table,
        rownames="gene_id"
    )

    write_tsv(de_res, "${id}_de_res_${condition_col}_edger.tsv")
    """
}


// process DE_DREAM {
//     /**
//      * For a pseodubulk analysis with mixed effects
//      * See https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/dream.html

//      * This is useful if there is a batch effect that is not controlled for
//      * by the experimental design
//      */

//     input:
//     tuple val(id), path(counts), path(samplesheet)
//     val(condition_col)
//     val(covariate_formula)

//     output:
//     path("de_res.csv"), emit: de_res

//     script:
//     """
//     # Run limma/dream as descirbed in their vignette
//     """

// }


process DE_MAST_MIXED_EFFECTS {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    cpus 40
    conda "/data/scratch/sturm/conda/envs/2020-pircher-mast"
    /**
     * Mixed effects model with MAST as described by Zimmermann et al.
     */

    input:
    tuple val(id), path(sce)
    val(condition_col)
    val(covariate_formula)

    output:
    path("*.tsv"), emit: de_res

    script:
    """
    #!/usr/bin/env Rscript

    options(mc.cores=${task.cpus})
    RhpcBLASctl::blas_set_num_threads(1)
    library(MAST)
    library(readr)

    sce = readRDS("${sce}")
    # expect counts in adata.X
    assays(sce)\$counts = assays(sce)\$X
    sce = scater::logNormCounts(sce)
    sca = MAST::SceToSingleCellAssay(sce, check_sanity = TRUE)

    res_zlm = zlm(
        # ~ leiden + n_genes_by_counts + (1 | organoid),
        ~ ${condition_col} + ${covariate_formula},
        sca,
        method="glmer",
        ebayes=FALSE,
        strictConvergence=FALSE,
        parallel=TRUE
    )

    tmp_contrast = "${condition_col}control"
    zlm_summary = MAST::summary(res_zlm, doLRT=tmp_contrast)

    summary_dt <- zlm_summary\$datatable
    de_res <- merge(summary_dt[summary_dt\$contrast==tmp_contrast
                                  & summary_dt\$component=='logFC', c(1,7,5,6,8)],
                        summary_dt[summary_dt\$contrast==tmp_contrast
                                  & summary_dt\$component=='H', c(1,4)],
                        by = 'primerid')

    write_tsv(de_res, "${id}_de_res_${condition_col}_mast.tsv")
    """
}


#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SPLIT_ANNDATA } from "../scconversion/main.nf"


process MAKE_PSEUDOBULK {
    cpus 1
    container "${baseDir}/containers/seuratdisk.sif"

    input:
    tuple val(id), path(input_adata)
    // the column containing the biological replicate variable
    val(replicate_col)
    val(condition_col)
    // minimum number of cells for a sample to be kept and
    // whether or not to keep unpaired samples
    tuple(val(min_cells_per_sample), val(keep_unpaired_samples))

    output:
    tuple val(id), path("${id}_counts.csv"), path("${id}_samplesheet.csv"), emit: pseudobulk

    script:
    def keep_unpaired_samples = keep_unpaired_samples ? "True" : "False"
    """
    #!/usr/bin/env python

    import scanpy as sc
    import pandas as pd
    import numpy as np

    adata = sc.read_h5ad("${input_adata}")
    replicate_col = "${replicate_col}"
    condition_col = "${condition_col}"
    keep_unpaired_samples = ${keep_unpaired_samples}
    min_cells_per_sample = ${min_cells_per_sample}

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
    if not keep_unpaired_samples:
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

process DE_DESEQ2 {
    /**
     * For standard pseudobulk analysis
     */
    cpus 2
    container "${baseDir}/containers/deseq2.sif"

    input:
    tuple val(id), path(counts), path(samplesheet)
    val(comparison)
    val(condition_col)
    val(covariate_formula)

    output:
    path("*_DESeq2_result.tsv"), emit: de_res

    when: comparison != "sum2zero"

    script:
    """
    run_deseq2.R $counts $samplesheet \\
        --cond_col $condition_col \\
        --c1 "${comparison[0]}" \\
        --c2 "${comparison[1]}" \\
        --resDir "." \\
        --prefix $id \\
        --covariate_formula "$covariate_formula" \\
        --cpus $task.cpus > ${id}_deseq2.log
    """
}


workflow deseq2_analysis {
    take:
        id                  // manually specified, unique identifier for the comparison
        adata
        column_to_test      // column to test, e.g. "tumor_stage"
        comparison          // Tuple ["treatment", "reference"] or String "sum2zero"
        cell_type_column    // column in adata containing cell-type information
        pseudobulk_group_by // column to generate pseudobulk by
        pseudobulk_settings // Tuple [min_cells, keep_unpaired samples]
        min_samples         // only consider cell-types with this minimum number of samples
        covariate_formula   // covariate formula (e.g. " + patient + sex")

    main:
    SPLIT_ANNDATA(
        adata,
        cell_type_column
    )
    MAKE_PSEUDOBULK(
        SPLIT_ANNDATA.out.adata.flatten().map{it -> [it.baseName, it]},
        pseudobulk_group_by,
        column_to_test,
        pseudobulk_settings
    )
    DE_DESEQ2(
        // only consider cell-types with at least N case/control samples
        MAKE_PSEUDOBULK.out.pseudobulk.filter{
            id, counts, samplesheet -> samplesheet.text.count("\n") >= min_samples
        }.map{
            //override ID
            it -> ["${id}_${it[0]}", it[1], it[2]]
        },
        comparison,
        column_to_test,
        covariate_formula
    )

    emit:
    deseq2_result = DE_DESEQ2.out

  }


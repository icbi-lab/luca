#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SPLIT_ANNDATA } from "../scconversion/main.nf"



process PREPARE_ANNDATA {
    cpus 1
    container "${baseDir}/containers/seuratdisk.sif"

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
    val(condition_col)
    val(covariate_formula)

    output:
    path("*_DESeq2_result.tsv"), emit: de_res

    script:
    """
    run_deseq2.R $counts $samplesheet \\
        --cond_col $condition_col \\
        --c1 case \\
        --c2 control \\
        --resDir "." \\
        --prefix $id \\
        --covariate_formula "$covariate_formula" \\
        --cpus $task.cpus > ${id}_deseq2.log
    """
}


workflow deseq2_analysis {
    take:
        adata
        column_to_test      // column to test, e.g. "tumor_stage"
        comparison          // Tuple [["treatment"], ["reference"]] or [["treatment"], "rest"]
        cell_type_column    // column in adata containing cell-type information
        pseudobulk_group_by // column to generate pseudobulkk by
        pseudobulk_settings // Tuple [min_cells, keep_unpaired samples]
        min_samples         // only consider cell-types with this minimum number of samples
        covariate_formula   // covariate formula (e.g. " + patient + sex")

    main:
    PREPARE_ANNDATA(
        adata,
        "X",
        column_to_test,
        comparison
    )
    SPLIT_ANNDATA(
        PREPARE_ANNDATA.out.adata,
        cell_type_column
    )
    MAKE_PSEUDOBULK(
        SPLIT_ANNDATA.out.adata.flatten().map{it -> [it.baseName, it]},
        pseudobulk_group_by,
        column_to_test,
        pseudobulk_settings
    )
    DE_DESEQ2(
        // only consider cell-types with at least three case/control samples
        MAKE_PSEUDOBULK.out.pseudobulk.filter{
            id, counts, samplesheet -> samplesheet.text.count("\n") >= min_samples
        },
        column_to_test,
        covariate_formula
    )

    emit:
    deseq2_result = DE_DESEQ2.out

  }


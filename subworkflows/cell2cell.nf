include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { JUPYTERNOTEBOOK as MAKE_CPDB_H5AD } from "../modules/local/jupyternotebook/main.nf"


process SQUIDPY {
    publishDir = [ path: { "${params.squidpy_outDir}"}, mode: 'link' ]

    input:
    tuple val(id), val(in_file)
    val(cell_type_key)

    output:
    path "*.pkl", emit: out_file, optional: true
    path "*.log", emit: log

    script:
    // using this instead of `errorStrategy` in order to also cache failed processes
    // (they will always fail due to characteristics of the data, e.g. too few cells)
    ignore_exit_code = task.ext.ignore_error ? "|| true" : ""
    """
    squidpy_cpdb.py \\
       -i ${in_file} \\
       -o ./ \\
       -c ${cell_type_key} \\
       -n ${task.cpus} \\
            > ${id}.log 2>&1 $ignore_exit_code
    """
}


workflow cell2cell {
    take: adata_annotated

    main:
    ch_adata_annotated = adata_annotated.map{ it -> [it.baseName, it]}
    SPLIT_ANNDATA(ch_adata_annotated, "sample")
    SQUIDPY(SPLIT_ANNDATA.out.adata.flatten().map{it -> [it.baseName, it]}, "cell_type_major")
    MAKE_CPDB_H5AD(
        Channel.value([[id: "61_make_cpdb_h5ad"], file("${baseDir}/analyses/60_cell2cell/61_make_cpdb_h5ad.py")]),
        [
            "adata_atlas": "full_atlas_merged.h5ad",
            "squidpy_dir": ".",
        ],
        adata_annotated.mix(SQUIDPY.out.out_file).collect()
    )

    emit:
    adata_cpdb = MAKE_CPDB_H5AD.out.artifacts
}

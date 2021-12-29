include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"


process SQUIDPY {
    publishDir = [ path: { "${params.squidpy_outDir}"}, mode: 'link' ]

    input:
    path in_file

    output:
    path "*", emit: out_file

    script:
    """
    squidpy_cpdb.py -i ${in_file} -o ./
    """
}


workflow cell2cell {
    take: adata_annotated

    main:
    ch_adata_annotated = Channel.value([adata_annotated.baseName, adata_annotated])
    SPLIT_ANNDATA(ch_adata_annotated, "sample")
    SQUIDPY(SPLIT_ANNDATA.out.adata.flatten())
}

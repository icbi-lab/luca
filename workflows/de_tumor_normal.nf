include { JUPYTERNOTEBOOK as PREPARE_FOR_DE }  from "../modules/local/jupyternotebook/main.nf"
include { SPLIT_ANNDATA as SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { PREPARE_ANNDATA as PREPARE_ANNDATA_TUMOR_NORMAL }  from "../modules/local/scde/main.nf"
include { MAKE_PSEUDOBULK as MAKE_PSEUDOBULK_TUMOR_NORMAL }  from "../modules/local/scde/main.nf"
include { DE_EDGER as DE_EDGER_TUMOR_NORMAL }  from "../modules/local/scde/main.nf"


workflow de_tumor_normal {
    take:
        adata_annotated

    main:
    PREPARE_FOR_DE(
        Channel.value([
            [id: "prepare_for_de"],
            file("${baseDir}/analyses/40_de_tumor_normal/41_prepare_de_analysis.py")
        ]),
        [
            "input_adata": "full_atlas_annotated.h5ad",
        ],
        adata_annotated
    )
    ch_adata_tumor_normal = PREPARE_FOR_DE.out.artifacts.map { it -> [it.baseName, it] }


    PREPARE_ANNDATA_TUMOR_NORMAL(
        ch_adata_tumor_normal,
        "X",
        "origin",
        [["tumor_primary"], "rest"]
    )
    SPLIT_ANNDATA(
        PREPARE_ANNDATA_TUMOR_NORMAL.out.adata,
        "cell_type"
    )
    MAKE_PSEUDOBULK_TUMOR_NORMAL(
        SPLIT_ANNDATA.out.adata.flatten().map{it -> [it.baseName, it]},
        "patient",
        "origin",
        [10, false]
    )
    DE_EDGER_TUMOR_NORMAL(
        // only consider cell-types with at least three case/control samples
        MAKE_PSEUDOBULK_TUMOR_NORMAL.out.pseudobulk.filter{
            id, counts, samplesheet -> samplesheet.text.count("\n") >= 6
        },
        "origin",
        "+ patient"
    )

}

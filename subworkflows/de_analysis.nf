include { JUPYTERNOTEBOOK as PREPARE_FOR_DE }  from "../modules/local/jupyternotebook/main.nf"
include {
    SPLIT_ANNDATA as SPLIT_ANNDATA_TUMOR_NORMAL;
    SPLIT_ANNDATA as SPLIT_ANNDATA_LUAD_LUSC;
    SPLIT_ANNDATA as SPLIT_ANNDATA_EARLY_ADVANCED }  from "../modules/local/scconversion/main.nf"
include {
    PREPARE_ANNDATA as PREPARE_ANNDATA_TUMOR_NORMAL;
    PREPARE_ANNDATA as PREPARE_ANNDATA_LUAD_LUSC;
    PREPARE_ANNDATA as PREPARE_ANNDATA_EARLY_ADVANCED }  from "../modules/local/scde/main.nf"
include {
    MAKE_PSEUDOBULK as MAKE_PSEUDOBULK_TUMOR_NORMAL;
    MAKE_PSEUDOBULK as MAKE_PSEUDOBULK_LUAD_LUSC ;
    MAKE_PSEUDOBULK as MAKE_PSEUDOBULK_EARLY_ADVANCED }  from "../modules/local/scde/main.nf"
include {
    DE_DESEQ2 as DE_DESEQ2_TUMOR_NORMAL;
    DE_DESEQ2 as DE_DESEQ2_EARLY_ADVANCED;
    DE_DESEQ2 as DE_DESEQ2_LUAD_LUSC }  from "../modules/local/scde/main.nf"


workflow de_analysis {
    take:
        adata_annotated

    main:
    PREPARE_FOR_DE(
        Channel.value([
            [id: "prepare_for_de"],
            file("${baseDir}/analyses/40_de_analysis/41_prepare_de_analysis.py")
        ]),
        adata_annotated.map{it -> [
            "input_adata": it.name,
        ]},
        adata_annotated
    )
    ch_prepare_adata = PREPARE_FOR_DE.out.artifacts.flatten().map { it -> [it.baseName, it] }

    /** TUMOR vs NORMAL ADJACENT samples (paired analysis) **/
    PREPARE_ANNDATA_TUMOR_NORMAL(
        ch_prepare_adata.filter{ id, adata -> id == "adata_tumor_normal"},
        "X",
        "origin",
        [["tumor_primary"], "rest"]
    )
    SPLIT_ANNDATA_TUMOR_NORMAL(
        PREPARE_ANNDATA_TUMOR_NORMAL.out.adata,
        "cell_type_major"
    )
    MAKE_PSEUDOBULK_TUMOR_NORMAL(
        SPLIT_ANNDATA_TUMOR_NORMAL.out.adata.flatten().map{it -> [it.baseName, it]},
        "patient",
        "origin",
        [10, false]
    )
    DE_DESEQ2_TUMOR_NORMAL(
        // only consider cell-types with at least three case/control samples
        MAKE_PSEUDOBULK_TUMOR_NORMAL.out.pseudobulk.filter{
            id, counts, samplesheet -> samplesheet.text.count("\n") >= 6
        },
        "origin",
        "+ patient"
    )

    /** LUAD vs. LUSC comparison of pirmary tumor samples **/
    PREPARE_ANNDATA_LUAD_LUSC(
        ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"},
        "X",
        "condition",
        [["LUAD"], ["LUSC"]]
    )
    SPLIT_ANNDATA_LUAD_LUSC(
        PREPARE_ANNDATA_LUAD_LUSC.out.adata,
        "cell_type_major"
    )
    MAKE_PSEUDOBULK_LUAD_LUSC(
        SPLIT_ANNDATA_LUAD_LUSC.out.adata.flatten().map{it -> [it.baseName, it]},
        "patient",
        "condition",
        [10, true]
    )

    DE_DESEQ2_LUAD_LUSC(
        // only consider cell-types with at least 10 samples
        // some samples will fail anyway due to not having a full rank matrix
        MAKE_PSEUDOBULK_LUAD_LUSC.out.pseudobulk.filter{
            id, counts, samplesheet -> samplesheet.text.count("\n") >= 10
        },
        "condition",
        "+ dataset"
    )

    /** Early vs. advanced comparison of pirmary tumor samples **/
    PREPARE_ANNDATA_EARLY_ADVANCED(
        ch_prepare_adata.filter{ id, adata -> id == "adata_primary_tumor"},
        "X",
        "tumor_stage",
        [["early"], ["advanced"]]
    )
    SPLIT_ANNDATA_EARLY_ADVANCED(
        PREPARE_ANNDATA_EARLY_ADVANCED.out.adata,
        "cell_type_major"
    )
    MAKE_PSEUDOBULK_EARLY_ADVANCED (
        SPLIT_ANNDATA_EARLY_ADVANCED.out.adata.flatten().map{it -> [it.baseName, it]},
        "patient",
        "tumor_stage",
        [10, true]
    )

    DE_DESEQ2_EARLY_ADVANCED(
        // only consider cell-types with at least 10 samples
        // some samples will fail anyway due to not having a full rank matrix
        MAKE_PSEUDOBULK_EARLY_ADVANCED.out.pseudobulk.filter{
            id, counts, samplesheet -> samplesheet.text.count("\n") >= 10
        },
        "tumor_stage",
        "+ dataset"
    )

}

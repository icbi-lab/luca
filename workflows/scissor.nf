include { SPLIT_ANNDATA }  from "../modules/local/scconversion/main.nf"
include { H5AD_TO_SCE }  from "../modules/local/scconversion/main.nf"
include { RMARKDOWNNOTEBOOK as SCISSOR } from "../modules/local/rmarkdownnotebook/main.nf"

workflow scissor {
    take: adata_integrated

    main:
    ch_adata_integrated = adata_integrated.map{ it -> [it.baseName, it]}
    SPLIT_ANNDATA(ch_adata_integrated, "patient")
    H5AD_TO_SCE(SPLIT_ANNDATA.out.adata.flatten().map{ it -> [it.baseName, it]})

    ch_sce = H5AD_TO_SCE.out.sce
    SCISSOR(
        Channel.value([
            [id: 'scissor'],
            file("${baseDir}/analyses/50_scissor/51_scissor_single_sample.Rmd")
        ]),
        ch_sce.map{
            id, sce -> [
                'id': id,
                 'sce_path': sce.name,
                 'survival_path': 'mmc1.xlsx',
                 'kras_mutation': 'kras_mutated.tsv',
                 'egfr_mutation': 'egfr_mutated.tsv',
                 'braf_mutation': 'braf_mutated.tsv',
                 'tcga_tpm_path': 'tcga-lung-primary.rds'
            ]
        },
        ch_sce.map{
            id, sce -> [
                sce,
                file("$baseDir/tables/tcga/mmc1.xlsx", checkIfExists: true),
                file("$baseDir/tables/tcga/kras_mutated.tsv", checkIfExists: true),
                file("$baseDir/tables/tcga/braf_mutated.tsv", checkIfExists: true),
                file("$baseDir/tables/tcga/egfr_mutated.tsv", checkIfExists: true),
                file("$baseDir/data/13_tcga/tcga-lung-primary.rds")
            ]
        }
    )
}

nextflow.enable.dsl = 2

include { nxfvars } from "./nxfvars.nf"

def testnxfvars() {
    return "echo test\n"
}

process P11_MERGE_ALL {
    publishDir "${options.publish_dir}", mode: "${params.publish_dir_mode}"

    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"

    input:
        path(samplesheet)
        path(adatas)

    output:
        path("merged_all.h5ad"), emit: adata
        path("*.html"), emit: notebook

    script:
    testnxfvars() <<
    """
    nxfvars execute ${baseDir}/analyses/20_integrate_scrnaseq_data/21_merge_all.py merge_all_report.html
    """
}

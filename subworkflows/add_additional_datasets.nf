include { check_samplesheet } from '../modules/local/check_samplesheet'
include { SCQC } from "../modules/local/scqc/main"
include { SCQC_MERGE_STATS } from "../modules/local/scqc_merge_stats/main.nf"

/**
 * Project new data onto the dataset using scANVI
 */
workflow  add_additional_datasets {

    main:
    ch_samples = Channel.from(check_samplesheet(params.additional_input, baseDir))

    SCQC(
        [
            file("${baseDir}/modules/local/scqc/scqc-notebook.py", checkIfExists: true),
            file("${baseDir}/modules/local/scqc/qc_plots.py", checkIfExists: true)
        ],
        ch_samples
    )
    SCQC_MERGE_STATS(SCQC.out.qc_stats.collect())
}

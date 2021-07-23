@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

def check_samplesheet(input_file) {
    def samples = []
    for(line in parseCsv(new FileReader(input_file))) {
        adata = file(line["input_adata"])
        assert adata.exists() : "File does not exist: ${line['input_adata']}"
        meta = line.toMap()
        meta.remove('input_adata')
        samples << [meta, adata]
    }
    return samples
}

/*
 * Reformat input samplesheet and check validity
 */
process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path samplesheet

    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/bactanti/bin/
    //TODO meaningful check
    """
    cp $samplesheet samplesheet.valid.csv
    # check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to get list of [ sample, single_end?, [ fastq_1, fastq_2 ] ]
def check_samplesheet_paths(LinkedHashMap row) {
    def meta = row.findAll { key, val -> key != "input_adata" }

    File adata_file = new File(row.input_adata)
    assert file.exists() : "File not found"

    def array = [ meta, file(row.input_adata, checkIfExists: true)]
    return array
}

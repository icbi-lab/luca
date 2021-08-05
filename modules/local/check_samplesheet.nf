@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

def check_samplesheet(input_file, base_dir) {
    def samples = []
    for(line in parseCsv(new FileReader(input_file))) {
        adata_path = "${base_dir}/${line['input_adata']}"
        adata = file(adata_path)
        assert adata.exists() : "File does not exist: ${adata_path}"
        meta = line.toMap()
        meta.remove('input_adata')
        samples << [meta, adata]
    }
    return samples
}

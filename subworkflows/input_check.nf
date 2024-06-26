//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv
        allele_number
        mode // "default", "merge", "nanopore"

    main:
        SAMPLESHEET_CHECK ( samplesheet, allele_number )
            .csv
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it, allele_number, mode) }
            .set { reads }

    emit:
        reads                                     // channel: [ val(meta), [ reads ] ]
        versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row, allele_number, mode) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.start_allele_1 = row.start_allele_1
    meta.end_allele_1 = row.end_allele_1
    if (allele_number == 2) {
      meta.start_allele_2 = row.start_allele_2
      meta.end_allele_2 = row.end_allele_2
    }
    if (mode == "nanopore") {
        meta.single_end = true
    } else if (mode == "default" || mode == "merge") {
        meta.single_end = false
    }
    
    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

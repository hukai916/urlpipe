//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    allele_number

    main:
    SAMPLESHEET_CHECK ( samplesheet, allele_number )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it, allele_number) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row, allele_number) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.length_cutoff_1_low = row.length_cutoff_1_low
    meta.length_cutoff_1_high = row.length_cutoff_1_high
    if (allele_number == 2) {
      meta.length_cutoff_2_low = row.length_cutoff_2_low
      meta.length_cutoff_2_high = row.length_cutoff_2_high
    }
    meta.single_end = row.single_end.toBoolean()

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

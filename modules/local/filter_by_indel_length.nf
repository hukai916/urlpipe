process FILTER_BY_INDEL_LENGTH {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/urlpipe:1.1"

    input:
    tuple val(meta), path(reads_input)
    path parse_cigar

    output:
    tuple val(meta), path(reads_input),   emit: reads
    path "stat.csv",                      emit: stat
    path  "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def indel_length_cutoff = task.ext.indel_length_cutoff ?: 10

    """
    touch stat.csv

    # step1: filter fastq reads by indel length



    # step2: provide BAM files for filtered fastq reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

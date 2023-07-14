process PARSE_CIGAR {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/urlpipe:1.0"

    input:
    tuple val(meta), path(reads)
    path bam

    output:
    tuple val(meta), path(reads), emit: reads_input
    path "parse_cigar_*.csv",     emit: parse_cigar                
    path  "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    scutls bam -lpiref --input $bam -o parse_cigar_${bam}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

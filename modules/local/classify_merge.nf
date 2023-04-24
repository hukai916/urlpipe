process CLASSIFY_MERGE {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/bbmap:0.1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("merged/*.fastq.gz"),     emit: reads
    tuple val(meta), path("non_merged/*.fastq.gz"), emit: reads_others
    path "stat/*.csv",                              emit: csv
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p merged non_merged stat

    bbmerge.sh in1=${prefix}_1.fastq.gz in2=${prefix}_2.fastq.gz out=merged/${prefix}.fastq.gz outu1=non_merged/${prefix}_1.fastq.gz outu2=non_merged/${prefix}_2.fastq.gz $args

    read_count_merge=\$(echo \$(zcat merged/${prefix}.fastq.gz | wc -l)/4|bc)
    read_count_non_merge=\$(echo \$(zcat non_merged/${prefix}_1.fastq.gz | wc -l)/4|bc)

    echo "${prefix},\$read_count_merge,\$read_count_non_merge" > stat/${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

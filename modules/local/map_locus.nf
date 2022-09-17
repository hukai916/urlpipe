process MAP_LOCUS {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("1a_map_locus/locus/*.fastq.gz"),     emit: reads_locus
    tuple val(meta), path("1a_map_locus/misprimed/*.fastq.gz"), emit: reads_misprimed
    tuple val(meta), path("1a_map_locus/problem/*.fastq.gz"),   emit: reads_problem
    path "map_locus_stat.tsv",                                  emit: stat
    path  "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 1a_map_locus/locus 1a_map_locus/misprimed 1a_map_locus/problem

    map_locus.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz 1a_map_locus/locus 1a_map_locus/misprimed 1a_map_locus/problem $args
    gzip 1a_map_locus/locus/*
    gzip 1a_map_locus/misprimed/*
    gzip 1a_map_locus/problem/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

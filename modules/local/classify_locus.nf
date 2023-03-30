process CLASSIFY_LOCUS {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("on_target_locus/*.fastq.gz"),     emit: reads_locus
    tuple val(meta), path("misprimed/*.fastq.gz"), emit: reads_misprimed
    tuple val(meta), path("off_target_locus/*.fastq.gz"),   emit: reads_problem
    path "stat/*.csv",                             emit: stat
    path  "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir on_target_locus misprimed off_target_locus stat
    touch stat/${prefix}.csv

    classify_locus.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz on_target_locus misprimed off_target_locus stat/${prefix}.csv $args
    gzip on_target_locus/*
    gzip misprimed/*
    gzip off_target_locus/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

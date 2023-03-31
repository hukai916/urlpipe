process CLASSIFY_LOCUS {
    // on_target_locus: both ends match
    // off_target_locus: neither ends match
    // problem_locus: one end matches
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("on_target_locus/*.fastq.gz"),  emit: reads_locus
    tuple val(meta), path("off_target_locus/*.fastq.gz"), emit: reads_problem
    tuple val(meta), path("problem_reads/*.fastq.gz"),    emit: reads_misprimed
    path "stat/*.csv",                                    emit: stat
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    classify_locus.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz on_target_locus off_target_locus problem_reads stat/${prefix}.csv $args

    gzip */*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

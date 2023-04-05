process REPEAT_LENGTH_DISTRIBUTION_DEFAULT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("raw_*.csv"),                 emit: raw_repeat_length_per_read_default
    tuple val(meta), path("repeat_length_per_*.csv"),   emit: repeat_length_per_read_default
    tuple val(meta), path("repeat_length_count_*.csv"), emit: repeat_length_count_default
    path "repeat_length_count_*.csv",                   emit: repeat_length_count_default_pure
    tuple val(meta), path("diagnosis_*.csv"),           emit: diagnosis_repeat_length_count_default
    path  "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"

    """
    # 1. output repeat_length_per_read_default_xxx.csv:
    repeat_length_per_read_default.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz raw_repeat_length_per_read_default_${prefix}.csv $args

    # 2. output repeat_length_count_xxx.csv and diagnosis_repeat_length_count_xxx.csv:
    repeat_length_count_default.py raw_repeat_length_per_read_default_${prefix}.csv repeat_length_per_read_default_${prefix}.csv repeat_length_count_${prefix}.csv diagnosis_repeat_length_count_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

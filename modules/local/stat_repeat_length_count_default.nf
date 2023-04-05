process STAT_REPEAT_LENGTH_COUNT_DEFAULT {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path stat

    output:
    path "*.csv",         emit: stat
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    stat_repeat_length_count_default.py *.csv repeat_length_count_default_umi_0.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

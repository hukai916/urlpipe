process STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    path csv

    output:
    path "*.csv",         emit: csv
    path "*.html",        emit: html
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    stat_repeat_length_count_default.py *.csv $args 500
    # 500 is the xlimit
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

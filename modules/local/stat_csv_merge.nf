process STAT_CSV_MERGE {
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
    stat_csv_merge.py *.csv "sample_id,read_count_merge,read_count_non_merge" classify_merge.csv classify_merge.html "" ".csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

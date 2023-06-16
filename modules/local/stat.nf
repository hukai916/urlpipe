process STAT {
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
    def outfile = task.ext.outfile ?: ''
    def header = task.ext.header ?: ''

    """
    (echo -e "$header" && cat *.csv | sort -n) > ${outfile}.txt
    mv ${outfile}.txt ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

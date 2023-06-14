process STAT_ADAPTOR {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(reads)
    // path reads
    path csv

    output:
    path "adaptor_counts.csv", emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    count_total=\$(expr \$(zcat $reads | wc -l) / 4)

    cat *.csv | awk -v total=\$count_total 'BEGIN{FS=","; OFS="," } { print \$1,\$2,\$2/total }' > adaptor_counts.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

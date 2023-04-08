process STAT_READ_COUNT_PER_UMI_CUTOFF {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    path csv
    val umi_cutoffs

    output:
    path "*.csv",         emit: csv
    path "*.html",        emit: html
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      stat_read_count_per_umi_cutoff.py read_count_*_umi_cutoff_\$i.csv read_count_*_umi_cutoff_\$i.html  "read_count_" "_umi_cutoff_\$i.csv"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

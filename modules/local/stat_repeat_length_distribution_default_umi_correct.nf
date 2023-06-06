process STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT {
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
      stat_repeat_length_count_default.py repeat_length_count_default_*_umi_\$i.csv repeat_length_count_default_umi_\$i.csv repeat_length_count_default_umi_\$i.html "repeat_length_count_default_" "_umi_\$i.csv" 500
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

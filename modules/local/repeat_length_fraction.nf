process REPEAT_LENGTH_FRACTION {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(csv_token)
    path csv_umi_0
    path csv_umi_x
    val umi_cutoffs

    output:
    path "*.csv",         emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      touch test.csv
      echo "test"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

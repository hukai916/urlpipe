process PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(csv)
    val umi_cutoffs

    output:
    path "*.html",        emit: html
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      plot_repeat_length_distribution_per_umi.py repeat_length_distribution_per_umi_\${i}_*.csv \$i repeat_length_distribution_per_umi_\${i}_${prefix}.html
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

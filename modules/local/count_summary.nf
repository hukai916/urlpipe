process COUNT_SUMMARY {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path frac_count_cutoff_0
    path frac_count_cutoff_x
    path indel_count_cutoff_x
    val umi_cutoffs
    val mode
    val outdir

    output:
    path "*/**/*.csv",     emit: count_summary
    path  "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/cutoff_\$i/${mode}
      count_summary.py *_frac_cutoff_\$i.csv *_indel_cutoff_$\i.csv ${outdir}/cutoff_\$i/${mode}/count_summary.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

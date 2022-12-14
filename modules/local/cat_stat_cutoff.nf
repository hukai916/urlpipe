process CAT_STAT_CUTOFF {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path stat
    val mode
    val header
    val umi_cutoffs
     // "mode", "mean", "ld"
    val output_prefix
    val outdir

    output:
    path "*/**/*_cutoff_*.csv",     emit: stat
    path  "versions.yml",           emit: versions

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
      (echo -e "$header" && cat *cutoff_\$i.csv | sort -n) > ${outdir}/cutoff_\$i/${mode}/${output_prefix}_cutoff_\$i.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

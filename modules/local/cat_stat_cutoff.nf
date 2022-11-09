process CAT_STAT_CUTOFF {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path stat
    val outdir
    val header
    val umi_cutoffs
    val mode // "mode", "mean", "ld"

    output:
    path "**/all_sample_cutoff_*.csv",     emit: stat
    path  "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/frac_\$i
      (echo -e "$header" && cat *cutoff_\$i.csv | sort -n) > ${outdir}/frac_\$i/${mode}/all_sample_cutoff_\$i.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

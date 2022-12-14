process PLOT_FRAC_CUTOFF {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path csv
    path stat_csv
    val umi_cutoffs
    val mode // mode, mean, ld
    val csv_prefix
    val outdir

    output:
    path "**/${mode}/*/all_sample_*.png",      emit: plot
    path  "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/${mode}/plot_frac_barplot
      mkdir -p ${outdir}/${mode}/plot_read_length_violin
      #csv: all_sample_cutoff_xxx.csv
      plot_frac.py ${csv_prefix}_cutoff_\$i.csv ${outdir}/${mode}/plot_frac_barplot/all_sample_frac_barplot_cutoff_\$i.png ${outdir}/${mode}/plot_read_length_violin/all_sample_read_length_violin_cutoff_\$i.png "\$i"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

process PLOT_FRAC_CUTOFF {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"



    PLOT_FRAC_CUTOFF_R1 (
      CAT_STAT5.out.stat, // all_sample_frac_cutoff_0.csv
      CAT_STAT_CUTOFF.out.stat, // all_sample_frac_cutoff_x.csv
      REPEAT_DIST_UMI_CORRECT_R1.out.stat_raw.collect(), // individual stat without UMI correction
      REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_mode_stat.collect(),
      // params.umi_cutoffs,
      params.allele_number,
      "mode",
      "all_sample_frac",
      "5d_r1_repeat_dist_umi_correct/plot_all_samples" // plot_read_length_violin and plot_frac_barplot
    )


    input:
    path all_sample_frac_cutoff_0_csv
    path all_sample_frac_cutoff_x_csv
    path stat_csv
    path stat_csv_cutoff
    val umi_cutoffs
    val csv_prefix
    val outdir

    output:
    path "**/${mode}/*/all_sample_*.png",      emit: plot
    path  "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    umi_cutoffs_str="0,$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/${mode}/plot_frac_barplot
      mkdir -p ${outdir}/${mode}/plot_read_length_violin

      plot_frac.py ${csv_prefix}_cutoff_\$i.csv ${outdir}/${mode}/plot_frac_barplot/all_sample_frac_barplot_cutoff_\$i.png ${outdir}/${mode}/plot_read_length_violin/all_sample_read_length_violin_cutoff_\$i.png "\$i"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

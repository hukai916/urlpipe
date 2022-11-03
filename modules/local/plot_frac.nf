process PLOT_FRAC {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path csv
    path stat_csv // raw stat csv, used for plot_frac.py though not mentioned in nf file
    val umi_cutoffs
    val outdir

    output:
    path "*/*/all_sample_*.png",  emit: plot
    path  "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    mkdir -p ${outdir}

    plot_frac.py $csv ${outdir}/all_sample_frac_barplot.png ${outdir}/all_sample_read_length_violin.png "$umi_cutoffs"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

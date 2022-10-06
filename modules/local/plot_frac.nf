process PLOT_FRAC {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path csv
    val outdir

    output:
    path "*/*/all_sample.png",        emit: plot
    path  "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    mkdir -p ${outdir}

    plot_frac.py $csv ${outdir}/all_sample.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

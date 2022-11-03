process PLOT_FRAC_UMI_CUTOFF {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path csv_raw // without correction
    tuple val(meta), path(csv_cutoff)
    val umi_cutoffs
    val outdir

    output:
    path "${outdir}/*.png", emit: plot
    path  "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    umi_cutoffs_str="$umi_cutoffs"

    plot_frac_umi_cutoff.py $csv_raw "$umi_cutoffs" $outdir/${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

process PLOT_UMI_GROUPS {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv_no_correction) // this is without correction
    path stat_csv // with corretion at different cutoffs
    val umi_cutoffs
    val outdir

    output:
    path "*/*/all_sample_*.png",          emit: plot
    path  "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}
      touch ${outdir}/all_sample_.npg
      #plot_umi_groups.py $csv_no_correction "$umi_cutoffs" ${outdir}/${prefix}.png
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

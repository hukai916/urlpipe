process REPEAT_DIST_WITHIN_UMI_GROUP {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(tsv)
    val outdir

    output:
    path "*/UMI_3/stat/*.tsv",        emit: umi_3_stat
    path "*/UMI_3/plot/*.png",        emit: umi_3_plot
    path "*/UMI_10/stat/*.tsv",       emit: umi_10_stat
    path "*/UMI_10/plot/*.png",       emit: umi_10_plot
    path "*/UMI_100/stat/*.tsv",      emit: umi_100_stat
    path "*/UMI_100/plot/*.png",      emit: umi_100_plot
    path  "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # examine repeat length distribution for UMI group 3, 10, and 100, each select 10 UMI to plot

    mkdir -p ${outdir}/UMI_3/stat ${outdir}/UMI_3/plot ${outdir}/UMI_10/stat ${outdir}/UMI_10/plot ${outdir}/UMI_100/stat ${outdir}/UMI_100/plot

    repeat_dist_within_umi_group.py $tsv $prefix 3 ${outdir}/UMI_3/stat ${outdir}/UMI_3/plot $args
    repeat_dist_within_umi_group.py $tsv $prefix 10 ${outdir}/UMI_10/stat ${outdir}/UMI_10/plot $args
    repeat_dist_within_umi_group.py $tsv $prefix 100 ${outdir}/UMI_100/stat ${outdir}/UMI_100/plot $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

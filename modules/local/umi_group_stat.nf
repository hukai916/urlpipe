process UMI_GROUP_STAT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    path stat
    val outdir

    output:
    tuple val(meta), path("${outdir}/*.csv"),       emit: stat
    path "*/input",                           emit: stat_raw
    path  "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # examine repeat length distribution for UMI group 3, 10, and 100, each select 10 UMI to plot

    mkdir -p ${outdir}/input
    cp $stat ${outdir}/input/ # otherwise, the emit will also place it under results/ directly

    umi_group_stat.py $csv $prefix ${outdir} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

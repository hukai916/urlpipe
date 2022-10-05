process REPEAT_DIST_DISTANCE_MERGED {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    path "*/stat/*.csv",                    emit: stat
    path "*/plot/*.png",                    emit: plot
    tuple val(meta), path("*/count/*.csv"), emit: count
    path "*/frac/*.csv",                    emit: frac
    path  "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${outdir}/stat ${outdir}/plot ${outdir}/count ${outdir}/frac

    repeat_dist_distance_merged.py ${prefix}.fastq.gz ${prefix} ${outdir} $args

    # calculate fractions
    tem="$args_frac"
    suffix="\${tem// /_}"

    calculate_frac.py $prefix ${outdir}/stat/${prefix}.csv ${outdir}/frac "$args_frac"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

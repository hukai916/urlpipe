process REPEAT_DIST_DISTANCE {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    path "*/stat_r1",              emit: stat_raw
    path "*/stat_r2",              emit: stat_r2
    path "*/plot_r1/*.png",        emit: plot_r1
    path "*/plot_r2/*.png",        emit: plot_r2
    tuple val(meta), path("*/count_r1/*.csv"),       emit: count_r1
    tuple val(meta), path("*/count_r2/*.csv"),       emit: count_r2
    path "*/frac_r1/*.csv",        emit: frac_r1
    path "*/frac_r2/*.csv",        emit: frac_r2
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${outdir}/stat_r1 ${outdir}/stat_r2 ${outdir}/plot_r1 ${outdir}/plot_r2 ${outdir}/count_r1 ${outdir}/count_r2

    repeat_dist_distance.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} ${outdir} $args

    # calculate fractions
    tem="$args_frac"
    suffix="\${tem// /_}"
    mkdir ${outdir}/frac_r1 ${outdir}/frac_r2

    calculate_frac.py $prefix ${outdir}/stat_r1/${prefix}.csv ${outdir}/frac_r1 "$args_frac"
    calculate_frac.py $prefix ${outdir}/stat_r2/${prefix}.csv ${outdir}/frac_r2 "$args_frac"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

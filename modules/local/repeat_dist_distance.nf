process REPEAT_DIST_DISTANCE {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    path "4d_repeat_distribution_distance/stat_r1/*.csv",        emit: stat_r1
    path "4d_repeat_distribution_distance/stat_r2/*.csv",        emit: stat_r2
    path "4d_repeat_distribution_distance/plot_r1/*.png",        emit: plot_r1
    path "4d_repeat_distribution_distance/plot_r2/*.png",        emit: plot_r2
    tuple val(meta), path("4d_repeat_distribution_distance/count_r1/*.csv"),       emit: count_r1
    tuple val(meta), path("4d_repeat_distribution_distance/count_r2/*.csv"),       emit: count_r2
    tuple val(meta), path("4d_repeat_distribution_distance/frac_r1_*/*.csv"),        emit: frac_r1
    tuple val(meta), path("4d_repeat_distribution_distance/frac_r2_*/*.csv"),        emit: frac_r2
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 4d_repeat_distribution_distance/stat_r1 4d_repeat_distribution_distance/stat_r2 4d_repeat_distribution_distance/plot_r1 4d_repeat_distribution_distance/plot_r2 4d_repeat_distribution_distance/count_r1 4d_repeat_distribution_distance/count_r2

    repeat_dist_distance.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} 4d_repeat_distribution_distance $args

    # calculate fractions
    suffix="\${$args_frac// /_}"
    mkdir 4d_repeat_distribution_distance/frac_r1 4d_repeat_distribution_distance/frac_r2

    calculate_frac.py $prefix 4d_repeat_distribution_distance/stat_r1/${prefix}.csv 4d_repeat_distribution_distance/frac_r1_\$suffix "$args_frac"
    calculate_frac.py $prefix 4d_repeat_distribution_distance/stat_r2/${prefix}.csv 4d_repeat_distribution_distance/frac_r2_\$suffix "$args_frac"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

process REPEAT_DIST_WITHIN_UMI_GROUP {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    path "4d_repeat_distribution_distance/stat_r1/*.tsv",        emit: stat_r1
    path "4d_repeat_distribution_distance/stat_r2/*.tsv",        emit: stat_r2
    path "4d_repeat_distribution_distance/plot_r1/*.png",        emit: plot_r1
    path "4d_repeat_distribution_distance/plot_r2/*.png",        emit: plot_r2
    path "4d_repeat_distribution_distance/count_r1/*.tsv",       emit: count_r1
    path "4d_repeat_distribution_distance/count_r2/*.tsv",       emit: count_r2
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 4d_repeat_distribution_distance/stat_r1 4d_repeat_distribution_distance/stat_r2 4d_repeat_distribution_distance/plot_r1 4d_repeat_distribution_distance/plot_r2 4d_repeat_distribution_distance/count_r1 4d_repeat_distribution_distance/count_r2

    repeat_dist_count.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} 4d_repeat_distribution_distance $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

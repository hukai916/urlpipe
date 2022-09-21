process REPEAT_DIST_DISTANCE_MERGED {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    path "4d_repeat_distribution_distance_merged/stat/*.csv",        emit: stat
    path "4d_repeat_distribution_distance_merged/plot/*.png",        emit: plot
    path "4d_repeat_distribution_distance_merged/count/*.csv",       emit: count
    path  "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 4d_repeat_distribution_distance_merged/stat 4d_repeat_distribution_distance_merged/plot 4d_repeat_distribution_distance_merged/count
    repeat_dist_distance_merged.py ${prefix}.fastq.gz ${prefix} 4d_repeat_distribution_distance_merged $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

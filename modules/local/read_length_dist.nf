process READ_LENGTH_DIST {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    path "*/stat_r1/*.stat.csv",   emit: stat_raw
    path "*/stat_r2/*.stat.csv",   emit: stat_r2
    path "*/plot_r1/*.png",        emit: plot_r1
    path "*/plot_r2/*.png",        emit: plot_r2
    tuple val(meta), path("*/count_r1/*.csv"),       emit: count_r1
    tuple val(meta), path("*/count_r2/*.csv"),       emit: count_r2
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    read_length_dist.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} ${outdir} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

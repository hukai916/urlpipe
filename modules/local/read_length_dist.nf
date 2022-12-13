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
    path "*/plot/*.png",        emit: plot
    tuple val(meta), path("*/count/*.csv"),       emit: count

    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${outdir}/stat_r1 ${outdir}/stat_r2 ${outdir}/plot_r1 ${outdir}/plot_r2 ${outdir}/count_r1 ${outdir}/count_r2

    read_length_dist.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} ${outdir} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

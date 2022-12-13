process READ_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    path reads
    val umi_cutoffs
    val outdir

    output:
    path "*/stat_r2/*.stat.csv",   emit: stat_r2
    path "*/plot_r1/*.png",        emit: plot_r1
    path "*/plot_r2/*.png",        emit: plot_r2
    tuple val(meta), path("*/count_r1/*.csv"),       emit: count_r1
    tuple val(meta), path("*/count_r2/*.csv"),       emit: count_r2
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    read_umi_correct.py $csv ${prefix} ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${outdir}/fastq/ "$umi_cutoffs"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

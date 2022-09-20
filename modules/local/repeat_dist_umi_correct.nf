process REPEAT_DIST_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(tsv)
    val outdir

    output:
    // path "*/cutoff_1/*.png",       emit: png_cutoff_1
    // path "*/cutoff_3/*.png",       emit: png_cutoff_3
    // path "*/cutoff_10/*.png",      emit: png_cutoff_10
    // path "*/cutoff_30/*.png",      emit: png_cutoff_30
    // path "*/cutoff_100/*.png",     emit: png_cutoff_100
    path "*/*", emit: test
    path  "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # plot repeat length distribution with UMI cutoff: 1, 3, 10, 30, 100

    mkdir -p ${outdir}/cutoff_1 ${outdir}/cutoff_3 ${outdir}/cutoff_10 ${outdir}/cutoff_30 ${outdir}/cutoff_100

    repeat_dist_umi_correct.py $tsv $prefix ${outdir}/cutoff_1 1
    repeat_dist_umi_correct.py $tsv $prefix ${outdir}/cutoff_3 3
    repeat_dist_umi_correct.py $tsv $prefix ${outdir}/cutoff_10 10
    repeat_dist_umi_correct.py $tsv $prefix ${outdir}/cutoff_30 30
    repeat_dist_umi_correct.py $tsv $prefix ${outdir}/cutoff_100 100

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

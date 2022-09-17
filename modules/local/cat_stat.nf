process CAT_STAT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path stat
    val header

    output:
    path "1a_map_locus/stat/all_sample.tsv",      emit: stat
    path  "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p 1a_map_locus/stat
    (echo $header && cat *.tsv) > 1a_map_locus/stat/all_sample.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

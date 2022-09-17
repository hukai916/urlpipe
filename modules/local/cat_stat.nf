process CAT_STAT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path stat
    val header

    output:
    path "1a_map_locus/stat/all.tsv",                             emit: stat
    path  "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 1a_map_locus/stat && touch 1a_map_locus/stat/all.tsv
    cat *.tsv > 1a_map_locus/stat/all_sample.tsv
    (echo $header && cat 1a_map_locus/stat/all_sample.tsv) > filename1 && mv filename1 1a_map_locus/stat/all_sample.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

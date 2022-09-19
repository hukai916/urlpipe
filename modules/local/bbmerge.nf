process BBMERGE {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/bbmap:0.1"

    input:
    tuple val(meta), path(reads)

    output:
    path "4b_bbmerge/merged/*.fastq.gz",        emit: reads_merged
    path "4b_bbmerge/non_merged/*.fastq.gz",    emit: reads_non_merged
    path "4b_bbmerge/stat/*.tsv",               emit: stat
    path  "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 4b_bbmerge/merged 4b_bbmerge/non_merged 4b_bbmerge/stat

    bbmerge.sh

    classify_readthrough.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz 4a_classify_readthrough/readthrough 4a_classify_readthrough/non_readthrough 4a_classify_readthrough/stat ${prefix} $args

    gzip 4a_classify_readthrough/readthrough*
    gzip 4a_classify_readthrough/non_readthrough/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

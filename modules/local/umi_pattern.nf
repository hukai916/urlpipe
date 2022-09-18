process UMI_PATTERN {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    path "2a_umi_pattern/stat/*.tsv",    emit: stat
    path "2a_umi_pattern/plot/*",        emit: plot
    path  "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 2a_umi_pattern/stat 2a_umi_pattern/plot

    zcat ${prefix}_1.fastq.gz | awk 'NR%4==1 {n = split(\$1, array, "_"); print array[n]}' | sort | uniq -c | awk '{print \$1}' | sort | uniq -c | awk '{print \$2, "\t", \$1}' | sort -n > tem.txt

    (echo -e "${prefix}\tcount" && cat tem.txt | sort -n) > 2a_umi_pattern/stat/${prefix}_UMI_distribution.tsv
    rm tem.txt

    plot_umi_pattern.py 2a_umi_pattern/stat/${prefix}_UMI_distribution.tsv ${prefix} 2a_umi_pattern/plot/${prefix}_UMI_distribution.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

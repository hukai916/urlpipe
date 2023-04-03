process READ_PER_UMI {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    path "*/stat/*.csv",    emit: stat
    path "*/plot/*",        emit: plot
    path  "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p stat plot

    zcat ${prefix}_1.fastq.gz | awk 'NR%4==1 {n = split(\$1, array, "_"); print array[n]}' | sort | uniq -c | awk '{print \$1}' | sort | uniq -c | awk '{print \$2, ",", \$1}' | sort -n > tem.txt

    (echo -e "${prefix},count" && cat tem.txt | sort -n) > stat/${prefix}_UMI_distribution.csv
    rm tem.txt

    plot_read_per_umi.py stat/${prefix}_read_per_umi.csv ${prefix} plot/${prefix}_read_per_umi.png $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

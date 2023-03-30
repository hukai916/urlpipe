process UMI_EXTRACT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/umitools_xenial:0.2"
    // ref: awk array: https://stackoverflow.com/questions/39703124/how-to-access-last-index-of-array-from-split-function-inside-awk
    // ref: let awk access bash variables: https://stackoverflow.com/questions/19075671/how-do-i-use-shell-variables-in-an-awk-script

    input:
    tuple val(meta), path(reads, stageAs: "input_fastq/*")

    output:
    tuple val(meta), path("*.fastq.gz"),  emit: reads
    path  "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    # extract UMI from start of the R1 read
    umi_tools extract $args -I input_fastq/${prefix}_1.fastq.gz -L log_${prefix}_1.txt -S ${prefix}_1.fastq.gz

    # match altered read name in R2, otherwise cutadapt complains
    zcat ${prefix}_1.fastq.gz > ${prefix}_1.fastq
    zcat input_fastq/${prefix}_2.fastq.gz > ${prefix}_2.fastq

    umi_length=\$(awk 'NR==1 {n = split(\$1, array, "_"); print length(array[n]); exit}' ${prefix}_1.fastq)

    awk -v umi_length="\$umi_length" 'NR == FNR { if (NR % 4 ==1) {umi_dict[substr(\$1, 1, length(\$1) - umi_length - 1)] = substr(\$1, length(\$1) - umi_length + 1, length(\$1))}; next } { if (NR %4 == 1) { print \$1 "_" umi_dict[\$1] " " \$2 } else { print \$0}}' ${prefix}_1.fastq ${prefix}_2.fastq | gzip > ${prefix}_2.fastq.gz

    rm *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UMI-tools: \$( umi_tools --version | sed -e "s/UMI-tools //g" )
    END_VERSIONS

    """
}

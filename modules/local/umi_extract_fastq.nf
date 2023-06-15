process UMI_EXTRACT_FASTQS {
    label 'process_low'

    container "hukai916/umi_tools:0.1"
    
    // ref: awk array: https://stackoverflow.com/questions/39703124/how-to-access-last-index-of-array-from-split-function-inside-awk
    // ref: let awk access bash variables: https://stackoverflow.com/questions/19075671/how-do-i-use-shell-variables-in-an-awk-script

    input:
    path reads

    output:
    path "umi_fastq_umi_to_keep/*.fastq.gz", emit: reads
    path  "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def umi_base_to_keep = task.ext.umi_base_to_keep ?: ''

    """
    mkdir umi_fastq umi_fastq_umi_to_keep
    for x in *.fastq.gz; do
        file_name=\${x%.fastq.gz}
        umi_tools extract $args -I \$x -L log_\${file_name}.txt -S umi_fastq/\$x

        umi_base_to_keep.py umi_fastq/\$x $umi_base_to_keep umi_fastq_umi_to_keep/\$x
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UMI-tools: \$( umi_tools --version | sed -e "s/UMI-tools //g" )
    END_VERSIONS

    """
}

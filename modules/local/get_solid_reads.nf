process GET_SOLID_READS {
    label 'process_low'

    container "hukai916/bwa:0.1"
    // ref: awk array: https://stackoverflow.com/questions/39703124/how-to-access-last-index-of-array-from-split-function-inside-awk
    // ref: let awk access bash variables: https://stackoverflow.com/questions/19075671/how-do-i-use-shell-variables-in-an-awk-script

    input:
    path reads
    path ref

    output:
    path "umi_fastq/*.fastq.gz",    emit: reads
    path  "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir 200bp_5p 200bp_3p bam_200bp_5p bam_200bp_3p
    bwa index $ref 

    for x in *.fastq.gz; do
        file_name=\${x%.fastq.gz}

        # Step1: cut off 200bp_5p and 200bp_3p
        zcat \$x | awk '{if (NR % 4 == 0) {print substr(\$0, 1, 200)} else if (NR % 4 == 2) {print substr(\$0, 1, 200)} else {print \$0} }' | gzip > 200bp_5p/\$x;
        zcat \$x | awk '{if (NR % 4 == 0) {print substr(\$0, length(\$0) - 199) } else if (NR % 4 == 2) {print  substr(\$0, length(\$0) - 199) } else {print \$0} }' | gzip > 200bp_3p

        # Step2: map against ref.fa using bwa
        bwa mem $ref 200bp_5p\$x | samtools view -bS -o bam_200bp_5p/\${file_name}.bam


        umi_tools extract $args -I \$x -L log_\${file_name}.txt -S umi_fastq/\$x
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UMI-tools: \$( umi_tools --version | sed -e "s/UMI-tools //g" )
    END_VERSIONS

    """
}

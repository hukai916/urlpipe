process BWA {
    label 'process_medium'
    container "hukai916/miniconda3_bwa:0.1"

    input:
    path ref 
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), path(reads), emit: bam_reads
    path  "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    def args = task.ext.args ?: ''

    """
    bwa index $ref
    bwa mem $args -t $task.cpus $ref $reads | samtools view -F 256 | samtools sort -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

process CLASSIFY_INDEL_NANOPORE {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/miniconda3_bio:1.0"

    input:
    tuple val(meta), path(reads)
    path ref

    output:
    tuple val(meta), path("no_indel/*.fastq.gz"),         emit: reads_no_indel
    tuple val(meta), path("indel_5p/*.fastq.gz"),         emit: reads_indel_5p
    tuple val(meta), path("indel_3p/*.fastq.gz"),         emit: reads_indel_3p
    tuple val(meta), path("indel_5p_and_3p/*.fastq.gz"),  emit: reads_indel_5p_3p
    tuple val(meta), path("indel_5p_or_3p/*.fastq.gz"),   emit: reads_indel_5p_or_3p
    path "indel_5p_or_3p/*.fastq.gz",                     emit: reads_indel_5p_or_3p_pure
    path "stat/*.csv",                                    emit: stat
    path "no_indel/*.bam",                                emit: bam_no_indel
    path "indel_5p/*.bam",                                emit: bam_indel_5p
    
    path "no_indel/*.bai",                                emit: bam_index_no_indel
    path  "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def indel_cutoff = task.ext.indel_cutoff ?: 0.5
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p no_indel indel_5p indel_3p indel_5p_3p stat

    classify_indel_nanopore.py ${prefix}.fastq.gz no_indel indel_5p indel_3p indel_5p_and_3p indel_5p_or_3p stat ${prefix} $args $indel_cutoff
    gzip */*.fastq

    # Add bam files and indices:
    bwa index $ref
    bwa mem -t $task.cpus $ref no_indel/*.fastq.gz | samtools view -bS | samtools sort -o no_indel/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_5p/*.fastq.gz | samtools view -bS | samtools sort -o indel_5p/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_3p/*.fastq.gz | samtools view -bS | samtools sort -o indel_3p/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_5p_and_3p/*.fastq.gz | samtools view -bS | samtools sort -o indel_5p_and_3p/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_5p_or_3p/*.fastq.gz | samtools view -bS | samtools sort -o indel_5p_or_3p/${prefix}.bam
    
    samtools index no_indel/*.bam
    samtools index indel_5p/*.bam
    samtools index indel_3p/*.bam
    samtools index indel_5p_and_3p/*.bam
    samtools index indel_5p_or_3p/*.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

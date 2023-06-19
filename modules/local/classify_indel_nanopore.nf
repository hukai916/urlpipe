process CLASSIFY_INDEL_NANOPORE {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/urlpipe:0.3"
    // container "hukai916/scutls:0.7"

    input:
    tuple val(meta), path(reads)
    path ref

    output:
    tuple val(meta), path("no_indel/*.fastq.gz"),         emit: reads_no_indel
    tuple val(meta), path("indel_5p_only/*.fastq.gz"),    emit: reads_indel_5p_only
    tuple val(meta), path("indel_3p_only/*.fastq.gz"),    emit: reads_indel_3p_only
    tuple val(meta), path("indel_5p_and_3p/*.fastq.gz"),  emit: reads_indel_5p_and_3p
    tuple val(meta), path("undetermined/*.fastq.gz"),     emit: reads_undetermined
    path "stat/*.csv",                                    emit: stat
    path "*/*.bam",                                       emit: bam
    path "*/*.bai",                                       emit: bam_index
    path  "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def repeat_flanking_left  = task.ext.repeat_flanking_left ?: ''
    def repeat_flanking_right = task.ext.repeat_flanking_right ?: ''
    def allowed_mismatch = task.ext.allowed_mismatch ?: ''
    

    // def args = task.ext.args ?: ''
    // def indel_cutoff = task.ext.indel_cutoff ?: 0.5
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p no_indel indel_5p_only indel_3p_only indel_5p_and_3p stat

    scutls barcode -l $repeat_flanking_left -nproc $task.cpus \\
        --input $reads \\
        -p 0 \\
        -e $allowed_mismatch
        -m > repeat_flanking_left.txt
    
    scutls barcode -l $repeat_flanking_right -nproc $task.cpus \\
        --input $reads \\
        -p -1 \\
        -e $allowed_mismatch
        -m > repeat_flanking_right.txt
    
    classify_indel_nanopore_location_info.py $reads $prefix \\
        repeat_flanking_left.txt repeat_flanking_right.txt \\
        no_indel/${prefix}.fastq.gz \\
        indel_5p_only/${prefix}.fastq.gz \\
        indel_3p_only/${prefix}.fastq.gz \\
        indel_5p_and_3p/${prefix}.fastq.gz \\
        undetermined/${prefix}.fastq.gz \\
        stat/${prefix}.csv

    # Add bam files and indices:
    bwa index $ref
    bwa mem -t $task.cpus $ref no_indel/*.fastq.gz | samtools view -bS | samtools sort -o no_indel/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_5p_only/*.fastq.gz | samtools view -bS | samtools sort -o indel_5p_only/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_3p_only/*.fastq.gz | samtools view -bS | samtools sort -o indel_3p_only/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_5p_and_3p/*.fastq.gz | samtools view -bS | samtools sort -o indel_5p_and_3p/${prefix}.bam
    bwa mem -t $task.cpus $ref undetermined/*.fastq.gz | samtools view -bS | samtools sort -o undetermined/${prefix}.bam

    samtools index no_indel/*.bam
    samtools index indel_5p_only/*.bam
    samtools index indel_3p_only/*.bam
    samtools index indel_5p_and_3p/*.bam
    samtools index undetermined/*.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

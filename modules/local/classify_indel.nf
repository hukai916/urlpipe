process CLASSIFY_INDEL {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/urlpipe:0.2"

    input:
    tuple val(meta), path(reads)
    path ref
    val ref_repeat_start
    val ref_repeat_end

    output:
    tuple val(meta), path("no_indel/*.fastq.gz"),         emit: reads_no_indel
    tuple val(meta), path("indel_5p/*.fastq.gz"),         emit: reads_indel_5p
    tuple val(meta), path("indel_3p/*.fastq.gz"),         emit: reads_indel_3p
    tuple val(meta), path("indel_5p_and_3p/*.fastq.gz"),  emit: reads_indel_5p_3p
    tuple val(meta), path("indel_5p_or_3p/*.fastq.gz"),   emit: reads_indel_5p_or_3p
    path "indel_5p_or_3p/*.fastq.gz",                     emit: reads_indel_5p_or_3p_pure
    path "stat/*.csv",                                    emit: stat
    path "*/*.bam",                                       emit: bam
    path "*/*.bai",                                       emit: bam_index
    path  "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ref_before_repeat_bp_to_check = task.ext.ref_before_repeat_bp_to_check ?: 20
    def ref_after_repeat_bp_to_check = task.ext.ref_after_repeat_bp_to_check ?: 20
    def m = task.ext.m ?: 1
    def indel_cutoff = task.ext.indel_cutoff ?: 0.5
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p no_indel indel_5p indel_3p indel_5p_3p stat

    start_tem=\$(($ref_repeat_start - 1))
    R1_tem=\$(get_seq.py $ref start \$start_tem no)

    end_tem=\$(($ref_repeat_end + $ref_after_repeat_bp_to_check))
    R2_tem=\$(get_seq.py $ref start \$end_tem no)
    echo ">ref1\n"\$R1_tem > ref1.fa
    echo ">ref2\n"\$R2_tem > ref2.fa

    R1=\$(get_seq.py ref1.fa end $ref_before_repeat_bp_to_check no)
    R2=\$(get_seq.py ref2.fa start $ref_after_repeat_bp_to_check no)

    classify_indel.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz no_indel indel_5p indel_3p indel_5p_and_3p indel_5p_or_3p stat ${prefix} \$R1 \$R2 $m $indel_cutoff

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

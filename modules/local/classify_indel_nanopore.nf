process CLASSIFY_INDEL_NANOPORE {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/urlpipe:0.4"
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
    path "stat/*.csv",                                    emit: stat_indel
    path "stat/mean_qc/mean_qc*",                         emit: mean_qc
    path "stat/per_site_qc/per_site_qc*.csv",             emit: per_site_qc
    path "*/bwa/*.bam",                                   emit: bam_bwa
    path "*/bwa/*.bai",                                   emit: bam_index_bwa
    path "indel_5p_only/minimap2/*.bam",                  emit: bam_minimap2_indel_5p_only
    path "indel_3p_only/minimap2/*.bam",                  emit: bam_minimap2_indel_3p_only
    path "indel_5p_and_3p/minimap2/*.bam",                emit: bam_minimap2_indel_5p_and_3p
    path "*/minimap2/*.bai",                              emit: bam_index_minimap2
    
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
    mkdir -p no_indel/bwa indel_5p_only/bwa indel_3p_only/bwa indel_5p_and_3p/bwa undetermined/bwa no_indel/minimap2 indel_5p_only/minimap2 indel_3p_only/minimap2 indel_5p_and_3p/minimap2 undetermined/minimap2 stat/mean_qc stat/per_site_qc

    scutls barcode -l $repeat_flanking_left -nproc $task.cpus \\
        --input $reads \\
        -p 0 \\
        -e $allowed_mismatch \\
        -m > repeat_flanking_left.txt
    
    scutls barcode -l $repeat_flanking_right -nproc $task.cpus \\
        --input $reads \\
        -p -1 \\
        -e $allowed_mismatch \\
        -m > repeat_flanking_right.txt
    
    classify_indel_nanopore_location_info.py $reads $prefix \\
        repeat_flanking_left.txt repeat_flanking_right.txt \\
        no_indel/${prefix}.fastq.gz \\
        indel_5p_only/${prefix}.fastq.gz \\
        indel_3p_only/${prefix}.fastq.gz \\
        indel_5p_and_3p/${prefix}.fastq.gz \\
        undetermined/${prefix}.fastq.gz \\
        stat/${prefix}.csv

    # Add bam files and indices using bwa:
    bwa index $ref
    bwa mem -t $task.cpus $ref no_indel/*.fastq.gz | samtools view -F 2048 -bS | samtools sort -o no_indel/bwa/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_5p_only/*.fastq.gz | samtools view -F 2048 -bS | samtools sort -o indel_5p_only/bwa/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_3p_only/*.fastq.gz | samtools view -F 2048 -bS | samtools sort -o indel_3p_only/bwa/${prefix}.bam
    bwa mem -t $task.cpus $ref indel_5p_and_3p/*.fastq.gz | samtools view -F 2048 -bS | samtools sort -o indel_5p_and_3p/bwa/${prefix}.bam
    bwa mem -t $task.cpus $ref undetermined/*.fastq.gz | samtools view -F 2048 -bS | samtools sort -o undetermined/bwa/${prefix}.bam

    samtools index no_indel/bwa/*.bam
    samtools index indel_5p_only/bwa/*.bam
    samtools index indel_3p_only/bwa/*.bam
    samtools index indel_5p_and_3p/bwa/*.bam
    samtools index undetermined/bwa/*.bam

    # Add bam files using minimap2
    minimap2 $ref no_indel/*.fastq.gz -a | samtools view -F 2048 -bS | samtools sort -o no_indel/minimap2/${prefix}.bam
    minimap2 $ref indel_5p_only/*.fastq.gz -a | samtools view -F 2048 -bS | samtools sort -o indel_5p_only/minimap2/${prefix}.bam
    minimap2 $ref indel_3p_only/*.fastq.gz -a | samtools view -F 2048 -bS | samtools sort -o indel_3p_only/minimap2/${prefix}.bam
    minimap2 $ref indel_5p_and_3p/*.fastq.gz -a | samtools view -F 2048 -bS | samtools sort -o indel_5p_and_3p/minimap2/${prefix}.bam
    minimap2 $ref undetermined/*.fastq.gz -a | samtools view -F 2048 -bS | samtools sort -o undetermined/minimap2/${prefix}.bam
    
    samtools index no_indel/minimap2/*.bam
    samtools index indel_5p_only/minimap2/*.bam
    samtools index indel_3p_only/minimap2/*.bam
    samtools index indel_5p_and_3p/minimap2/*.bam
    samtools index undetermined/minimap2/*.bam

    # Add read_quality plot1: mean overall read quality per read category
    get_mean_qc.py no_indel/${prefix}.fastq.gz ${prefix} stat/mean_qc/mean_qc_no_indel_${prefix}.csv
    get_mean_qc.py indel_5p_only/${prefix}.fastq.gz ${prefix} stat/mean_qc/mean_qc_indel_5p_only_${prefix}.csv
    get_mean_qc.py indel_3p_only/${prefix}.fastq.gz ${prefix} stat/mean_qc/mean_qc_indel_3p_only_${prefix}.csv
    get_mean_qc.py indel_5p_and_3p/${prefix}.fastq.gz ${prefix} stat/mean_qc/mean_qc_indel_5p_and_3p_${prefix}.csv
    get_mean_qc.py undetermined/${prefix}.fastq.gz ${prefix} stat/mean_qc/mean_qc_undetermined_${prefix}.csv

    # Add read quality plot2: mean quality per mapped site per read category
    samtools mpileup no_indel/minimap2/${prefix}.bam -o stat/per_site_qc/per_site_qc_no_indel_${prefix}.txt
    samtools mpileup indel_5p_only/minimap2/${prefix}.bam -o stat/per_site_qc/per_site_qc_indel_5p_only_${prefix}.txt
    samtools mpileup indel_3p_only/minimap2/${prefix}.bam -o stat/per_site_qc/per_site_qc_indel_3p_only_${prefix}.txt
    samtools mpileup indel_5p_and_3p/minimap2/${prefix}.bam -o stat/per_site_qc/per_site_qc_indel_5p_and_3p_${prefix}.txt
    samtools mpileup undetermined/minimap2/${prefix}.bam -o stat/per_site_qc/per_site_qc_undetermined_${prefix}.txt

    cat stat/per_site_qc/per_site_qc_no_indel_${prefix}.txt | awk 'BEGIN {OFS=","} {print "${prefix}",\$2,\$6}' > stat/per_site_qc/per_site_qc_no_indel_${prefix}.csv
    cat stat/per_site_qc/per_site_qc_indel_5p_only_${prefix}.txt | awk 'BEGIN {OFS=","} {print "${prefix}",\$2,\$6}' > stat/per_site_qc/per_site_qc_indel_5p_only_${prefix}.csv
    cat stat/per_site_qc/per_site_qc_indel_3p_only_${prefix}.txt | awk 'BEGIN {OFS=","} {print "${prefix}",\$2,\$6}' > stat/per_site_qc/per_site_qc_indel_3p_only_${prefix}.csv
    cat stat/per_site_qc/per_site_qc_indel_5p_and_3p_${prefix}.txt | awk 'BEGIN {OFS=","} {print "${prefix}",\$2,\$6}' > stat/per_site_qc/per_site_qc_indel_5p_and_3p_${prefix}.csv
    cat stat/per_site_qc/per_site_qc_undetermined_${prefix}.txt | awk 'BEGIN {OFS=","} {print "${prefix}",\$2,\$6}' > stat/per_site_qc/per_site_qc_undetermined_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

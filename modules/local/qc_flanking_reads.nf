process QC_FLANKING_READS {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/urlpipe:0.8"

    input:
    tuple val(meta), path(reads)
    path ref

    output:
    tuple val(meta), path(reads),                         emit: reads_input
    tuple val(meta), path("stat/read_id_mean_qc_*.csv"),  emit: read_id_mean_qc 
    path "stat/*",                                        emit: stat_all
    path "*/minimap2/*.bam",                              emit: bam_minimap2
    path "*/minimap2/*.bai",                              emit: bam_index_minimap2
    path  "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def left_flanking_coordinates  = task.ext.left_flanking_coordinates ?: '0:20'
    def right_flanking_coordinates = task.ext.right_flanking_coordinates ?: '100:120'
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir minimap2 stat

    # Step1: map using minimap2:
    minimap2 $ref $reads -a | samtools view -F 2048 -bS | samtools sort -o minimap2/${prefix}.bam
    samtools index minimap2/*.bam

    # Step2: parse out mean read quality for left_flanking_coordinates for each mapped read
    IFS=':' read -r left_s left_e <<< "$left_flanking_coordinates"
    scutls bam --input minimap2/${prefix}.bam -lpir \$left_s -o left_start.txt
    scutls bam --input minimap2/${prefix}.bam -lpir \$left_e -o left_end.txt
    get_fastq_mean_qc_by_range.py $reads left_start.txt left_end.txt left_flanking_qc.txt

    # Step3: parse out mean read quality for right_flanking_coordinates for each mapped read
    IFS=':' read -r right_s right_e <<< "$right_flanking_coordinates"
    scutls bam --input minimap2/${prefix}.bam -lpir \$right_s -o right_start.txt
    scutls bam --input minimap2/${prefix}.bam -lpir \$right_e -o right_end.txt
    get_fastq_mean_qc_by_range.py $reads right_start.txt right_end.txt right_flanking_qc.txt

    # Step4: plot mean quality distribution for left_flanking_qc.txt and right_flanking_qc.txt
    # also output read_id_mean_qc_*.csv
    plot_mean_qc_flanking.py left_flanking_qc.txt right_flanking_qc.txt stat/flanking_qc_${prefix}.html $prefix stat/read_id_mean_qc_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

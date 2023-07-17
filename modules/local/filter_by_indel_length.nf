process FILTER_BY_INDEL_LENGTH {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/urlpipe:1.1"

    input:
    tuple val(meta), path(reads_input)
    path parse_cigar

    output:
    tuple val(meta), path("indel_pass_filter/*.fastq.gz"), emit: reads_indel_pass_filter
    tuple val(meta), path("indel_not_pass_filter/*.fastq.gz"), emit: reads_indel_not_pass_filter
    path "stat/*_stat.csv",               emit: stat
    path "*/bwa/*.bam",                   emit: bam_bwa
    path "*/bwa/*.bai",                   emit: bam_index_bwa
    path "*/minimap2/*.bam",              emit: bam_minimap2
    path "*/minimap2/*.bai",              emit: bam_index_minimap2
    path  "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def indel_length_cutoff = task.ext.indel_length_cutoff ?: 10
    def ref_range = task.ext.ref_range ?: "1081:1101"

    """
    mkdir -p indel_pass_filter/minimap2 indel_not_pass_filter/minimap2 stat

    # step1: filter fastq reads by indel length
    get_fastq_by_parse_cigar.py $reads $parse_cigar $indel_length_cutoff $ref_range indel_pass_filter/$reads indel_not_pass_filter/$reads

    # step2: provive stats
    count_indel_pass_filter_reads=\$(get_fastq_count.py indel_pass_filter/$reads)
    count_indel_not_pass_filter_reads=\$(get_fastq_count.py indel_not_pass_filter/$reads)
    percent_indel_pass_filter_reads=\$(echo "scale=2; \$count_indel_pass_filter_reads / (\$count_indel_pass_filter_reads + \$count_indel_not_pass_filter_reads + 0.001)" | bc)
    
    percent_indel_not_pass_filter_reads=\$(echo "scale=2; \$count_indel_not_pass_filter_reads / (\$count_indel_pass_filter_reads + \$count_indel_not_pass_filter_reads + 0.001)" | bc)

    reads_tem=$reads
    filename=\${reads_tem%.fastq.gz}

    echo \$filename,\$count_indel_pass_filter_reads,\$count_indel_not_pass_filter_reads,\$percent_indel_pass_filter_reads,\$percent_indel_not_pass_filter_reads > stat/\${filename}_stat.csv

    # step3: add BAM files
    minimap2 $ref indel_pass_filter/$reads -a | samtools view -F 2048 -bS | samtools sort -o indel_pass_filter/minimap2/\${filename}.bam

    minimap2 $ref indel_not_pass_filter/$reads -a | samtools view -F 2048 -bS | samtools sort -o indel_not_pass_filter/minimap2/\${filename}.bam

    samtools index indel_pass_filter/minimap2/*.bam
    samtools index indel_not_pass_filter/minimap2/*.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

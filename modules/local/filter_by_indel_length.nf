process FILTER_BY_INDEL_LENGTH {
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/urlpipe:1.1"

    input:
    tuple val(meta), path(reads)
    path parse_cigar
    path ref

    output:
    tuple val(meta), path("indel_above_length_cutoff/*.fastq.gz"), emit: indel_above_length_cutoff
    tuple val(meta), path("indel_below_length_cutoff/*.fastq.gz"), emit: indel_below_length_cutoff
    path "stat/*_stat.csv",               emit: stat
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
    mkdir -p indel_above_length_cutoff/minimap2 indel_below_length_cutoff/minimap2 stat

    # step1: filter fastq reads by indel length
    get_fastq_seq_by_parse_cigar.py $reads $parse_cigar $indel_length_cutoff $ref_range indel_above_length_cutoff/$reads indel_below_length_cutoff/$reads

    # step2: calculate stats
    count_indel_above_length_cutoff=\$(get_fastq_count.py indel_above_length_cutoff/$reads)
    count_indel_below_length_cutoff=\$(get_fastq_count.py indel_below_length_cutoff/$reads)
    percent_indel_above_length_cutoff=\$(echo "scale=2; \$count_indel_above_length_cutoff / (\$count_indel_above_length_cutoff + \$count_indel_below_length_cutoff_reads + 0.001)" | bc)
    percent_indel_below_length_cutoff=\$(echo "scale=2; \$count_indel_below_length_cutoff / (\$count_indel_above_length_cutoff + \$count_indel_below_length_cutoff_reads + 0.001)" | bc)

    reads_tem=$reads
    filename=\${reads_tem%.fastq.gz}

    echo \$filename,\$count_indel_above_length_cutoff,\$count_indel_below_length_cutoff,\$percent_indel_above_length_cutoff_reads,\$percent_indel_below_length_cutoff_reads > stat/\${filename}_stat.csv

    # step3: add BAM files
    minimap2 $ref indel_above_length_cutoff/$reads -a | samtools view -F 2048 -bS | samtools sort -o indel_above_length_cutoff/minimap2/\${filename}.bam

    minimap2 $ref indel_below_length_cutoff/$reads -a | samtools view -F 2048 -bS | samtools sort -o indel_below_length_cutoff/minimap2/\${filename}.bam

    samtools index indel_above_length_cutoff/minimap2/*.bam
    samtools index indel_below_length_cutoff/minimap2/*.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

process GET_HIGH_QUALITY_FLANKING_READS {
    label 'process_medium'

    container "hukai916/urlpipe:0.9"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(read_id_mean_qc)

    output:
    path "high_quality/*.fastq.gz",           emit: reads
    path "low_quality/not_pass_*.fastq.gz",   emit: reads_low_quality_flanking
    path "stat/*_stat.csv",                   emit: stat

    path  "versions.yml",                     emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mean_quality_cutoff = task.ext.mean_quality_cutoff ?: '25'

    """
    mkdir -p high_quality low_quality stat

    # step1: filter reads based on mean quality score
    get_high_quality_flanking_reads.py $reads $mean_quality_cutoff mean_qc_cutoff high_quality/${reads} low_quality/${reads}

    # step2: obtain some statistics
    count_high_quality_reads=\$(get_fastq_count.py high_quality/$reads)
    count_low_quality_reads=\$(get_fastq_count.py low_quality/$reads)
    percent_high_quality_reads=\$(echo "scale=2; \$count_high_quality_reads / (\$count_high_quality_reads + \$count_low_quality_reads) + 0.01" | bc)
    percent_low_quality_reads=\$(echo "scale=2; \$count_low_quality_reads / (\$count_high_quality_reads + \$count_low_quality_reads) + 0.01" | bc)

    reads_tem=$reads
    filename=\${reads_tem%.fastq.gz}

    echo \$filename,\$count_high_quality_flanking_reads,\$count_low_quality_flanking_reads,\$percent_high_quality_flanking_reads,\$percent_low_quality_flanking_reads > stat/\${filename}_stat.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UMI-tools: \$( umi_tools --version | sed -e "s/UMI-tools //g" )
    END_VERSIONS

    """
}

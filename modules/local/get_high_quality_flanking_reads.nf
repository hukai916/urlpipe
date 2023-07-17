process GET_HIGH_QUALITY_FLANKING_READS {
    label 'process_medium'

    container "hukai916/urlpipe:1.0"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(read_id_mean_qc)
    tuple val(meta), path(read_id_mean_qc_extend)
    
    output:
    tuple val(meta), path("high_quality/*.fastq.gz"), emit: reads
    tuple val(meta), path("low_quality/*.fastq.gz"),  emit: reads_low_quality_flanking
    path "stat/*_stat.csv",                           emit: stat

    path  "versions.yml",                             emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mean_quality_cutoff = task.ext.mean_quality_cutoff ?: '25'
    def mode = task.ext.mode ?: 'default'

    """
    mkdir -p high_quality low_quality stat

    # step1: filter reads based on mean quality score
    get_high_quality_flanking_reads.py $reads $read_id_mean_qc $mean_quality_cutoff high_quality/${reads} low_quality/${reads} $mode $read_id_mean_qc_extend

    # step2: obtain some statistics
    count_high_quality_reads=\$(get_fastq_count.py high_quality/$reads)
    count_low_quality_reads=\$(get_fastq_count.py low_quality/$reads)
    percent_high_quality_reads=\$(echo "scale=2; \$count_high_quality_reads / (\$count_high_quality_reads + \$count_low_quality_reads + 0.001)" | bc)
    percent_low_quality_reads=\$(echo "scale=2; \$count_low_quality_reads / (\$count_high_quality_reads + \$count_low_quality_reads + 0.001)" | bc)

    reads_tem=$reads
    filename=\${reads_tem%.fastq.gz}

    echo \$filename,\$count_high_quality_reads,\$count_low_quality_reads,\$percent_high_quality_reads,\$percent_low_quality_reads > stat/\${filename}_stat.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

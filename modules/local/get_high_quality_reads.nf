process GET_HIGH_QUALITY_READS {
    label 'process_medium'

    container "hukai916/urlpipe:0.7"

    input:
    path reads

    output:
    path "high_quality/*.fastq.gz",           emit: reads
    path "low_quality/not_pass_*.fastq.gz",   emit: reads_low_quality
    path "stat/*_stat.csv",                   emit: stat

    path  "versions.yml",                     emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def quality_cutoff = task.ext.quality_cutoff ?: '25'

    """
    mkdir -p high_quality low_quality stat

    # step1: filter reads based on mean quality score
    scutls fastq --input $reads --output high_quality/$reads -fq $quality_cutoff

    # step2: obtain some statistics
    count_high_quality_reads=\$(expr \$(zcat high_quality/$reads | wc -l) / 4)
    count_low_quality_reads=\$(expr \$(zcat high_quality/not_pass_$reads | wc -l) / 4)
    percent_high_quality_reads=\$(echo "scale=2; \$count_high_quality_reads / (\$count_high_quality_reads + \$count_low_quality_reads)" | bc)
    percent_low_quality_reads=\$(echo "scale=2; \$count_low_quality_reads / (\$count_high_quality_reads + \$count_low_quality_reads)" | bc)

    reads_tem=$reads
    filename=\${reads_tem%.fastq.gz}

    echo \$filename,\$count_high_quality_reads,\$count_low_quality_reads,\$percent_high_quality_reads,\$percent_low_quality_reads > stat/\${filename}_stat.csv

    mv high_quality/not_pass_${reads} low_quality/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UMI-tools: \$( umi_tools --version | sed -e "s/UMI-tools //g" )
    END_VERSIONS

    """
}

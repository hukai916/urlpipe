process SPLIT_ALLELE {
    label 'process_cpu'

    container "hukai916/scutls:0.6"

    input:
    path reads

    output:
    path "allele_1/*.fastq.gz",     emit: reads_allele2
    path "allele_2/*.fastq.gz",     emit: reads_allele2
    path "undetermined/*.fastq.gz", emit: reads_undetermined
    path "stat/*_stat.csv",         emit: stat
    path "versions.yml",           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def snp_left_flanking = task.ext.snp_left_flanking ?: ''
    def snp_right_flanking = task.ext.snp_right_flanking ?: '' 
    def snp_1 = task.ext.snp_1 ?: ''
    def snp_2 = task.ext.snp_2 ?: ''
    def allowed_error = task.ext.allowed_error ?: '3'


    """
    mkdir full_length_reads partial_length_reads stat

    # step1: obtain bases from ref start and end
    ref_start=\$(get_fasta_range.py $ref $ref_start_range)
    ref_end=\$(get_fasta_range.py $ref $ref_end_range)

    # step2: obtain location in the read of ref_start and ref_start, as well as the read length
    scutls barcode -l \$ref_start -nproc $task.cpus \\
        --input $reads \\
        -p 0 \\
        -e $allowed_error > ref_start_in_range.txt

    scutls barcode -l \$ref_end -nproc $task.cpus \\
    --input $reads \\
    -p -1 \\
    -e $allowed_error > ref_end_in_range.txt
    
    # step3: obtain read length
    get_full_length_reads.py $reads ref_start_in_range.txt ref_end_in_range.txt $read_start_range $read_end_range full_length_reads/$reads partial_length_reads/$reads

    # step4: get some stats
    count_full_length_reads=\$(expr \$(zcat full_length_reads/*.fastq.gz | wc -l) / 4)
    count_partial_length_reads=\$(expr \$(zcat partial_length_reads/*.fastq.gz | wc -l) / 4)
    percent_full_length_reads=\$(echo "scale=2; \$count_full_length_reads / (\$count_full_length_reads + \$count_partial_length_reads)" | bc)
    percent_partial_length_reads=\$(echo "scale=2; \$count_partial_length_reads / (\$count_full_length_reads + \$count_partial_length_reads)" | bc)
    
    reads_tem=$reads
    filename=\${reads_tem%.fastq.gz}

    echo \$filename,\$count_full_length_reads,\$count_partial_length_reads,\$percent_full_length_reads,\$percent_partial_length_reads > stat/\${filename}_stat.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UMI-tools: \$( umi_tools --version | sed -e "s/UMI-tools //g" )
    END_VERSIONS

    """
}

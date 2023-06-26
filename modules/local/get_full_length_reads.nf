process GET_FULL_LENGTH_READS {
    label 'process_cpu'

    container "hukai916/urlpipe:0.6"

    input:
    path reads
    path ref

    output:
    path "full_length_reads/*.fastq.gz",    emit: reads
    path "partial_length_reads/*.fastq.gz", emit: reads_partial_length
    path "stat/*_stat.csv",                 emit: stat
    path "*/bwa/*.bam",                     emit: bam_bwa
    path "*/bwa/*.bai",                     emit: bam_index_bwa
    path "*/minimap2/*.bam",                emit: bam_minimap2
    path "*/minimap2/*.bai",                emit: bam_index_minimap2
    
    path  "versions.yml",                   emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ref_start_range = task.ext.ref_start_range ?: '0:50'
    def ref_end_range = task.ext.ref_end_range ?: '-50:' 
    def read_start_range = task.ext.read_start_range ?: '0:250'
    def read_end_range = task.ext.read_end_range ?: '-150:'
    def allowed_error = task.ext.allowed_error ?: '3'


    """
    mkdir -p full_length_reads/bwa full_length_reads/minimap2 partial_length_reads/bwa partial_length_reads/minimap2 stat

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

    # step5: add bam files
    ## Add bam files and indices using bwa:
    bwa index $ref
    bwa mem -t $task.cpus $ref full_length_reads/$reads | samtools view -F 2048 -bS | samtools sort -o full_length_reads/bwa/\${filename}.bam
    bwa mem -t $task.cpus $ref partial_length_reads/$reads | samtools view -F 2048 -bS | samtools sort -o partial_length_reads/bwa/\${filename}.bam

    samtools index full_length_reads/bwa/*.bam
    samtools index partial_length_reads/bwa/*.bam

    ## Add bam files using minimap2
    minimap2 $ref full_length_reads/$reads -a | samtools view -F 2048 -bS | samtools sort -o full_length_reads/minimap2/\${filename}.bam
    minimap2 $ref partial_length_reads/$reads -a | samtools view -F 2048 -bS | samtools sort -o partial_length_reads/minimap2/\${filename}.bam

    samtools index full_length_reads/minimap2/*.bam
    samtools index partial_length_reads/minimap2/*.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UMI-tools: \$( umi_tools --version | sed -e "s/UMI-tools //g" )
    END_VERSIONS

    """
}

process SPLIT_ALLELE {
    label 'process_cpu'

    container "hukai916/scutls:0.8"

    input:
    path reads

    output:
    path "allele_1/*.fastq.gz",     emit: reads_allele1
    path "allele_2/*.fastq.gz",     emit: reads_allele2
    path "undetermined/*.fastq.gz", emit: reads_undetermined
    path "stat/*_stat.csv",         emit: stat
    path "versions.yml",            emit: versions

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
    mkdir split_allele stat

    # step1: obtain SNP information from each read
    snp1="($snp_left_flanking){e<=$allowed_error}($snp_1)($snp_right_flanking){e<=$allowed_error}"
    snp2="($snp_left_flanking){e<=$allowed_error}($snp_2)($snp_right_flanking){e<=$allowed_error}"

    scutls barcode --locate "\$snp1" -i $reads -nproc $task.cpus > snp1.csv
    scutls barcode --locate "\$snp2" -i $reads -nproc $task.cpus > snp2.csv

    # step2: split reads into snp1, snp2, undetermined reads
    split_allele.py $reads snp1.csv snp2.csv split_allele/snp1_$reads split_allele/snp2_$reads split_allele/undetermined_$reads

    # step3: obtain some statistics
    snp1_reads=\$(expr \$(zcat split_allele/snp1_*.fastq.gz | wc -l) / 4)
    snp2_reads=\$(expr \$(zcat split_allele/snp2_*.fastq.gz | wc -l) / 4)
    undetermined_reads=\$(expr \$(zcat split_allele/undetermined_*.fastq.gz | wc -l) / 4)

    percent_snp1_reads=\$(echo "scale=2; \$snp1_reads / (\$snp1_reads + \$snp2_reads + \$undetermined_reads)" | bc)
    percent_snp2_reads=\$(echo "scale=2; \$snp2_reads / (\$snp1_reads + \$snp2_reads + \$undetermined_reads)" | bc)
    percent_undetermined_reads=\$(echo "scale=2; \$undetermined_reads / (\$snp1_reads + \$snp2_reads + \$undetermined_reads)" | bc)

    reads_tem=$reads
    filename=\${reads_tem%.fastq.gz}

    echo \$filename,\$snp1_reads,\$snp2_reads,\$undetermined_reads,\$percent_snp1_reads,\$percent_snp2_reads,\$percent_undetermined_reads > stat/\${filename}_stat.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UMI-tools: \$( umi_tools --version | sed -e "s/UMI-tools //g" )
    END_VERSIONS

    """
}

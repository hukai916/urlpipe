process GET_VALID_NANOPORE_READS {
    label 'process_cpu'

    container "hukai916/scutls:0.4"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_valid.fastq.gz"), emit: reads_valid 
    tuple val(meta), path("*_invalid.fastq.gz"), emit: reads_invalid
    tuple val(meta), path("*_valid_rc.fastq.gz"), emit: reads_valid_rc
    tuple val(meta), path("*_invalid_rc.fastq.gz"), emit: reads_invalid_rc
    tuple val(meta), path("*_valid_combine.fastq.gz"), emit: reads_valid_combine

    path "individual_csv/*.csv", emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # using original barcode:
    scutls barcode -c $args -nproc $task.cpus \\
        --input $reads \\
        -o ${prefix}_valid.fastq.gz \\
        -o2 ${prefix}_invalid.fastq.gz 

    # using rc of barcode:
    scutls barcode -rcb -c $args -nproc $task.cpus \\
    --input $reads \\
    -o ${prefix}_valid_rc.fastq.gz \\
    -o2 ${prefix}_invalid_rc.fastq.gz 

    # reverse complement the rc reads:
    seqtk seq -r ${prefix}_valid_rc.fastq.gz | gzip > ${prefix}_valid_rc_rc.fastq.gz

    # combine and remove potential redundancies
    zcat ${prefix}_valid.fastq.gz ${prefix}_valid_rc_rc.fastq.gz | gzip > tem.fastq.gz
    scutls fastq -i tem.fastq.gz -o ${prefix}_valid_combine.fastq.gz -u

    

    # some stats
    mkdir individual_csv
    count_valid=\$(expr \$(zcat ${prefix}_valid.fastq.gz | wc -l) / 4)
    count_invalid=\$(expr \$(zcat ${prefix}_invalid.fastq.gz | wc -l) / 4)
    count_valid_rc=\$(expr \$(zcat ${prefix}_valid_rc.fastq.gz | wc -l) / 4)
    count_invalid_rc=\$(expr \$(zcat ${prefix}_invalid_rc.fastq.gz | wc -l) / 4)
    count_valid_combine=\$(expr \$(zcat ${prefix}_valid_combine.fastq.gz | wc -l) / 4)
    
    echo ${prefix}_valid,\$count_valid > individual_csv/${prefix}_valid.csv
    echo ${prefix}_invalid,\$count_invalid > individual_csv/${prefix}_invalid.csv
    echo ${prefix}_valid_rc,\$count_valid_rc > individual_csv/${prefix}_valid_rc.csv
    echo ${prefix}_invalid_rc,\$count_invalid_rc > individual_csv/${prefix}_invalid_rc.csv
    echo ${prefix}_valid_combine,\$count_valid_combine > individual_csv/${prefix}_valid_combine.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

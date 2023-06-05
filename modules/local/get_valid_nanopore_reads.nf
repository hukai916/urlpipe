process GET_VALID_NANOPORE_READS {
    label 'process_cpu'

    container "hukai916/scutls:0.2"

    input:
    tuple val(meta), path(reads)

    output:
    path "*_valid.fastq.gz", emit: reads_valid 
    path "*_invalid.fastq.gz", emit: reads_invalid
    path "*_valid_rc.fastq.gz", emit: reads_valid_rc
    path "*_invalid_rc.fastq.gz", emit: reads_invalid_rc
    path "individual_csv/*.csv", emit: count_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when
    def prefix = task.ext.prefix ?: "${meta.id}"

    script:
    def args = task.ext.args ?: ''

    """
    # using original barcode:
    scutls barcode -c $args -nproc $task.cpus \\
        --input $reads \\
        -o ${prefix}_valid.fastq.gz \\
        -o2 ${prefix}_invalid.fastq.gz 

    # using rc of barcode
    scutls barcode -rcb -c $args -nproc $task.cpus \\
    --input $reads \\
    -o ${prefix}_valid_rc.fastq.gz \\
    -o2 ${prefix}_invalid_rc.fastq.gz 

    # some stats
    mkdir individual_csv
    count_valid=\$(expr \$(zcat ${prefix}_valid.fastq.gz | wc -l) / 4)
    count_invalid=\$(expr \$(zcat ${prefix}_invalid.fastq.gz | wc -l) / 4)
    count_valid_rc=\$(expr \$(zcat ${prefix}_valid_rc.fastq.gz | wc -l) / 4)
    count_invalid_rc=\$(expr \$(zcat ${prefix}_invalid_rc.fastq.gz | wc -l) / 4)
    echo ${prefix}_valid,\$count_valid > individual_csv/${prefix}_valid.csv
    echo ${prefix}_invalid,\$count_invalid > individual_csv/${prefix}_invalid.csv
    echo ${prefix}_valid_rc,\$count_valid_rc > individual_csv/${prefix}_valid_rc.csv
    echo ${prefix}_invalid_rc,\$count_invalid_rc > individual_csv/${prefix}_invalid_rc.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

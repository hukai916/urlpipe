process BARCODE_COUNT {
    label 'process_cpu'

    container "hukai916/scutls:0.1"

    input:
    tuple val(meta), path(reads)
    val prefix

    output:
    path "individual_fastq/*.fastq.gz", emit: fastq 
    path "individual_csv/*_with_bc.csv", emit: count_with_bc 
    path "individual_csv/*_without_bc.csv", emit: count_without_bc 
    path "individual_csv/*_with_bc_rc.csv", emit: count_with_bc_rc 
    path "individual_csv/*_without_bc_rc.csv", emit: count_without_bc_rc
    path "individual_csv/*.csv", emit: count_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir individual_fastq
    mkdir individual_csv

    # using original barcode:
    scutls barcode -c $args -nproc $task.cpus \\
        --input $reads \\
        -o individual_fastq/${prefix}_with_bc.fastq.gz \\
        -o2 individual_fastq/${prefix}_without_bc.fastq.gz 

    # using rc of barcode
    scutls barcode -rcb -c $args -nproc $task.cpus \\
    --input $reads \\
    -o individual_fastq/${prefix}_with_bc_rc.fastq.gz \\
    -o2 individual_fastq/${prefix}_without_bc_rc.fastq.gz 

    # some stats
    count_with_bc=\$(expr \$(zcat individual_fastq/${prefix}_with_bc.fastq.gz | wc -l) / 4)
    count_without_bc=\$(expr \$(zcat individual_fastq/${prefix}_without_bc.fastq.gz | wc -l) / 4)
    count_with_bc_rc=\$(expr \$(zcat individual_fastq/${prefix}_with_bc_rc.fastq.gz | wc -l) / 4)
    count_without_bc_rc=\$(expr \$(zcat individual_fastq/${prefix}_without_bc_rc.fastq.gz | wc -l) / 4)
    echo ${prefix}_with_bc,\$count_with_bc > individual_csv/${prefix}_with_bc.csv
    echo ${prefix}_without_bc,\$count_without_bc > individual_csv/${prefix}_without_bc.csv
    echo ${prefix}_with_bc_rc,\$count_with_bc_rc > individual_csv/${prefix}_with_bc_rc.csv
    echo ${prefix}_without_bc_rc,\$count_without_bc_rc > individual_csv/${prefix}_without_bc_rc.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

process BARCODE_COUNT {
    label 'process_cpu'

    container "hukai916/scutls:0.1"

    input:
    tuple val(meta), path(reads)
    val prefix

    output:
    path "*.fastq.gz", emit: fastq 
    path "*_with_bc.csv", emit: count_with_bc 
    path "*_without_bc.csv", emit: count_without_bc 
    path "*_with_bc_rc.csv", emit: count_with_bc_rc 
    path "*_without_bc_rc.csv", emit: count_without_bc_rc

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # using original barcode:
    scutls barcode -c $args -nproc $task.cpus \\
        --input $reads \\
        -o ${prefix}_with_bc.fastq.gz \\
        -o2 ${prefix}_without_bc.fastq.gz 

    # using rc of barcode
    scutls barcode -rcb -c $args -nproc $task.cpus \\
    --input $reads \\
    -o ${prefix}_with_bc_rc.fastq.gz \\
    -o2 ${prefix}_without_bc_rc.fastq.gz 

    # some stats
    count_with_bc=\$(expr \$(zcat ${prefix}_with_bc.fastq.gz | wc -l) / 4)
    count_without_bc=\$(expr \$(zcat ${prefix}_without_bc.fastq.gz | wc -l) / 4)
    count_with_bc_rc=\$(expr \$(zcat ${prefix}_with_bc_rc.fastq.gz | wc -l) / 4)
    count_without_bc_rc=\$(expr \$(zcat ${prefix}_without_bc_rc.fastq.gz | wc -l) / 4)
    echo ${prefix},\$count_with_bc > ${prefix}_with_bc.csv
    echo ${prefix},\$count_without_bc > ${prefix}_without_bc.csv
    echo ${prefix},\$count_with_bc_rc > ${prefix}_with_bc_rc.csv
    echo ${prefix},\$count_without_bc_rc > ${prefix}_without_bc_rc.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

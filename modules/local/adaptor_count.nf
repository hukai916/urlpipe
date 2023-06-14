process ADAPTOR_COUNT {
    label 'process_cpu'

    container "hukai916/scutls:0.3"

    input:
    tuple val(meta), path(reads)
    val prefix

    output:
    path "individual_fastq/*.fastq.gz", emit: fastq 
    path "individual_csv/*_with_ap.csv", emit: count_with_ap 
    path "individual_csv/*_without_ap.csv", emit: count_without_ap 
    path "individual_csv/*_with_ap_rc.csv", emit: count_with_ap_rc 
    path "individual_csv/*_without_ap_rc.csv", emit: count_without_ap_rc
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
        -o individual_fastq/${prefix}_with_ap.fastq.gz \\
        -o2 individual_fastq/${prefix}_without_ap.fastq.gz 

    # using rc of barcode
    scutls barcode -rcb -c $args -nproc $task.cpus \\
    --input $reads \\
    -o individual_fastq/${prefix}_with_ap_rc.fastq.gz \\
    -o2 individual_fastq/${prefix}_without_ap_rc.fastq.gz 

    # some stats
    count_with_ap=\$(expr \$(zcat individual_fastq/${prefix}_with_ap.fastq.gz | wc -l) / 4)
    count_without_ap=\$(expr \$(zcat individual_fastq/${prefix}_without_ap.fastq.gz | wc -l) / 4)
    count_with_ap_rc=\$(expr \$(zcat individual_fastq/${prefix}_with_ap_rc.fastq.gz | wc -l) / 4)
    count_without_ap_rc=\$(expr \$(zcat individual_fastq/${prefix}_without_ap_rc.fastq.gz | wc -l) / 4)
    echo ${prefix}_with_ap,\$count_with_ap > individual_csv/${prefix}_with_ap.csv
    echo ${prefix}_without_ap,\$count_without_ap > individual_csv/${prefix}_without_ap.csv
    echo ${prefix}_with_ap_rc,\$count_with_ap_rc > individual_csv/${prefix}_with_ap_rc.csv
    echo ${prefix}_without_ap_rc,\$count_without_ap_rc > individual_csv/${prefix}_without_ap_rc.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

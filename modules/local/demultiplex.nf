process DEMULTIPLEX {
    label 'process_medium'

    container "hukai916/pheniqs:0.2"

    input:
    tuple val(meta), path(reads)
    val reverse_complement
    val barcode_pos

    output:
    path "*.fastq.gz", emit: reads 
    path "*_pheniqs_report.json", emit: pheniqs_log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """

    # Step1: prepare Pheniqs json file
    prep_pheniqs_json.py $reads $prefix $reverse_complement $barcode_pos $args

    # Step2: perform Pheniqs demultiplexing
    pheniqs mux --config *.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

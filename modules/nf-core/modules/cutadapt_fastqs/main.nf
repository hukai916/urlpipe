process CUTADAPT_FASTQS {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? 'bioconda::cutadapt=3.4' : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
    //     'quay.io/biocontainers/cutadapt:3.4--py39h38f01e4_1' }"
    container "hukai916/cutadapt_xenial:0.1"

    input:
    path reads

    output:
    path 'trimmed/*.fastq.gz', emit: reads
    path '*.log',              emit: log
    path "versions.yml",       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir trimmed
    for x in *.fastq.gz; do
        cutadapt --cores $task.cpus \\
            $args -o trimmed/\$x \$x \\
            > \${x}.cutadapt.log
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}

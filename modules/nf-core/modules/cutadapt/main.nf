process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? 'bioconda::cutadapt=3.4' : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
    //     'quay.io/biocontainers/cutadapt:3.4--py39h38f01e4_1' }"
    container "hukai916/cutadapt_xenial:0.1"

    input:
    tuple val(meta), path(reads, stageAs: "input_fastq/*")

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log'),      emit: log
    path "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def trimmed  = meta.single_end ? "-o ${prefix}.fastq.gz" : "-o ${prefix}_1.fastq.gz -p ${prefix}_2.fastq.gz"
    """
    cutadapt \\
        --cores $task.cpus \\
        $trimmed \\
        $reads \\
        $args \\
        > ${prefix}.cutadapt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}

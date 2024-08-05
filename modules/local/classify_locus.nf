process CLASSIFY_LOCUS {
    // on_target_locus: both ends match
    // off_target_locus: neither ends match
    // problem_locus: one end matches
    tag "$meta.id"
    label 'process_medium'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    path ref

    output:
    tuple val(meta), path("on_target_locus/*.fastq.gz"),  emit: reads_locus
    tuple val(meta), path("off_target_locus/*.fastq.gz"), emit: reads_problem
    tuple val(meta), path("problem_reads/*.fastq.gz"),    emit: reads_misprimed
    path "stat/*.csv",                                    emit: stat
    path "versions.yml",                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ref_start_bp_to_check = task.ref_start_bp_to_check ?: 10
    def ref_end_bp_to_check = task.ref_end_bp_to_check ?: 10
    def m = task.m ?: 1
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    R1=\$(get_seq.py $ref start $ref_start_bp_to_check no)
    R2=\$(get_seq.py $ref end $ref_end_bp_to_check rc)

    classify_locus.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz on_target_locus off_target_locus problem_reads stat/${prefix}.csv \$R1 \$R2 $m

    gzip */*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

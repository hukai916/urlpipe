process REPEAT_LENGTH_DISTRIBUTION_DEFAULT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(reads)
    path ref
    val ref_repeat_start
    val ref_repeat_end

    output:
    tuple val(meta), path("raw_*.csv"),                 emit: raw_repeat_length_per_read_default
    tuple val(meta), path("repeat_length_per_*.csv"),   emit: repeat_length_per_read_default
    tuple val(meta), path("repeat_length_count_*.csv"), emit: repeat_length_count_default
    path "repeat_length_count_*.csv",                   emit: repeat_length_count_default_pure
    path "raw_repeat_length_per_read_*.csv",            emit: raw_repeat_length_per_read
    tuple val(meta), path("diagnosis_*.csv"),           emit: diagnosis_repeat_length_count_default
    path  "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ref_before_repeat_bp_to_check = task.ext.ref_before_repeat_bp_to_check ?: 20
    def ref_after_repeat_bp_to_check = task.ext.ref_after_repeat_bp_to_check ?: 20
    def m = task.ext.m ?: 1
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"

    """
    # 1. output repeat_length_per_read_default_xxx.csv:
    start_tem=\$(($ref_repeat_start - 1))
    R1_tem=\$(get_seq.py $ref start \$start_tem no)

    end_tem=\$(($ref_repeat_end + $ref_after_repeat_bp_to_check))
    R2_tem=\$(get_seq.py $ref start \$end_tem no)
    echo ">ref1\n"\$R1_tem > ref1.fa
    echo ">ref2\n"\$R2_tem > ref2.fa

    R1=\$(get_seq.py ref1.fa end $ref_before_repeat_bp_to_check no)
    R2=\$(get_seq.py ref2.fa end $ref_after_repeat_bp_to_check no)

    repeat_length_per_read_default.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz raw_repeat_length_per_read_default_${prefix}.csv \$R1 \$R1 $m

    # 2. output repeat_length_count_xxx.csv and diagnosis_repeat_length_count_xxx.csv:
    repeat_length_distribution_default.py raw_repeat_length_per_read_default_${prefix}.csv repeat_length_per_read_default_${prefix}.csv repeat_length_count_${prefix}.csv diagnosis_repeat_length_count_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

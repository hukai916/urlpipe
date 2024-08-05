process CLASSIFY_READTHROUGH {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    file ref 
    val ref_repeat_start
    val ref_repeat_end

    output:
    tuple val(meta), path("readthrough/*.fastq.gz"),     emit: reads_through
    tuple val(meta), path("non_readthrough/*.fastq.gz"), emit: reads_nonethrough
    path "stat/*.csv",    emit: stat
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ref_before_repeat_bp_to_check = task.ext.ref_before_repeat_bp_to_check ?: 20
    def ref_after_repeat_bp_to_check = task.ext.ref_after_repeat_bp_to_check ?: 20
    def m = task.ext.m ?: 1
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir readthrough non_readthrough stat

    start_tem=\$(($ref_repeat_start - 1))
    R1_tem=\$(get_seq.py $ref start \$start_tem no)

    end_tem=\$(($ref_repeat_end + $ref_after_repeat_bp_to_check))
    R2_tem=\$(get_seq.py $ref start \$end_tem no)
    echo ">ref1\n"\$R1_tem > ref1.fa
    echo ">ref2\n"\$R2_tem > ref2.fa

    R1=\$(get_seq.py ref1.fa end $ref_before_repeat_bp_to_check no)
    R2=\$(get_seq.py ref2.fa end $ref_after_repeat_bp_to_check no)

    classify_readthrough.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz readthrough non_readthrough stat ${prefix} \$R1 \$R2 $m

    gzip */*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

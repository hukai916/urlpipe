process BWA_LENGTH {
    label 'process_medium'
    container "hukai916/miniconda3_bio:0.4"

    input:
    tuple val(meta), path(bam), path(reads)

    output:
    tuple val(meta), path("raw_*.csv"),                 emit: raw_repeat_length_per_read_default
    tuple val(meta), path("repeat_length_per_*.csv"),   emit: repeat_length_per_read_default
    tuple val(meta), path("repeat_length_count_*.csv"), emit: repeat_length_count_default
    path "repeat_length_count_*.csv",                   emit: repeat_length_count_default_pure
    path "raw_repeat_length_per_read_*.csv",            emit: raw_repeat_length_per_read
    tuple val(meta), path("diagnosis_*.csv"),           emit: diagnosis_repeat_length_count_default
    path  "versions.yml",                               emit: versions

    script:
    def prefix = "${meta.id}"

    """
    # 1. output repeat_length_per_read_default_xxx.csv:
    samtools index $bam
        #repeat_length_per_read_bwa.py $bam $reads raw_repeat_length_per_read_bwa.csv
    repeat_length_per_read_bwa.py $bam $reads raw_repeat_length_per_read_default_${prefix}.csv
    
    # 2. output repeat_length_count_xxx.csv and diagnosis_repeat_length_count_xxx.csv:
        #repeat_length_distribution_default.py raw_repeat_length_per_read_default.csv repeat_length_per_read_default.csv repeat_length_count_default.csv diagnosis_repeat_length_count_default.csv
    repeat_length_distribution_default.py raw_repeat_length_per_read_default_${prefix}.csv repeat_length_per_read_default_${prefix}.csv repeat_length_count_${prefix}.csv diagnosis_repeat_length_count_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

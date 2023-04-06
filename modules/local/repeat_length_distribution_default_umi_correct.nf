process REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(csv)
    val umi_correction_method
    val umi_cutoffs

    output:
    tuple val(meta), path("UMI_ReadCount_ReadLengthCorrected_*"), emit: umi_readcount_readlength_corrected
    path "repeat_length_count_default_*", emit: repeat_length_count_default_umi_correct
    path  "versions.yml",                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # 1. obtain master table: UMI,read_count,repeat_length_corrected (by corresponding method)
    repeat_length_umi_correct.py $csv $umi_correction_method UMI_ReadCount_ReadLengthCorrected_${prefix}.csv

    # 2. obtain corrected table using various UMI cutoffs:
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      repeat_length_distribution_default_umi_correct.py UMI_ReadCount_ReadLengthCorrected_${prefix}.csv \$i repeat_length_count_default_${prefix}_umi_\$i.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

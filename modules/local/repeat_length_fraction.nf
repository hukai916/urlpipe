process REPEAT_LENGTH_FRACTION {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(csv_token)
    path csv_umi_0
    path csv_umi_x
    val allele_number
    val umi_cutoffs

    output:
    path "*.csv",         emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def start_allele_1  = "${meta.start_allele_1}"
    def end_allele_1 = "${meta.end_allele_1}"
    def start_allele_2  = "${meta.start_allele_2}" ?: 1 // simply place holder
    def end_allele_2 = "${meta.end_allele_2}" ?: 2

    """
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      repeat_length_frac.py ${allele_number} ${prefix} repeat_length_count_default_umi_\$i.csv repeat_length_fraction_umi_\$i.csv $start_allele_1 $end_allele_1 $start_allele_2 $end_allele_2
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

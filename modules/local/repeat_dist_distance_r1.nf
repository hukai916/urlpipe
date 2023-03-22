process REPEAT_DIST_DISTANCE_R1 {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val allele_number
    val outdir

    output:
    path "*/stat_r1/*.stat.csv",                emit: stat_raw
    path "*/plot_r1/*.png",                     emit: plot_r1
    tuple val(meta), path("*/count_r1/*.csv"),  emit: count_r1
    path "*/frac_r1/*.csv",                     emit: frac_r1
    path  "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def args_frac = task.ext.args_frac ?: ''
    def length_cutoff_1_low  = "${meta.length_cutoff_1_low}"
    def length_cutoff_1_high = "${meta.length_cutoff_1_high}"
    def length_cutoff_2_low  = "${meta.length_cutoff_2_low}" ?: 1 // simply place holder
    def length_cutoff_2_high = "${meta.length_cutoff_2_high}" ?: 2

    def prefix = "${meta.id}"

    """
    mkdir -p ${outdir}/stat_r1 ${outdir}/plot_r1 ${outdir}/count_r1 ${outdir}/frac_r1

    repeat_dist_distance_r1.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} ${outdir} $args

    # calculate fractions
    if [ $allele_number -eq 1 ]; then
      calculate_frac.py $prefix ${outdir}/stat_r1/${prefix}.stat.csv ${outdir}/frac_r1 0 $length_cutoff_1_low $length_cutoff_1_high
    elif [ $allele_number -eq 2 ]; then
      calculate_frac_2.py $prefix ${outdir}/stat_r1/${prefix}.stat.csv ${outdir}/frac_r1 0 $length_cutoff_1_low $length_cutoff_1_high $length_cutoff_2_low $length_cutoff_2_high
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

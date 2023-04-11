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
    def start_allele_1  = "${meta.start_allele_1}"
    def end_allele_1 = "${meta.end_allele_1}"
    def start_allele_2  = "${meta.start_allele_2}" ?: 1 // simply place holder
    def end_allele_2 = "${meta.end_allele_2}" ?: 2

    def prefix = "${meta.id}"

    """
    mkdir -p ${outdir}/stat_r1 ${outdir}/plot_r1 ${outdir}/count_r1 ${outdir}/frac_r1

    repeat_dist_distance_r1.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} ${outdir} $args

    # calculate fractions
    if [ $allele_number -eq 1 ]; then
      calculate_frac.py $prefix ${outdir}/stat_r1/${prefix}.stat.csv ${outdir}/frac_r1 0 $start_allele_1 $end_allele_1
    elif [ $allele_number -eq 2 ]; then
      calculate_frac_2.py $prefix ${outdir}/stat_r1/${prefix}.stat.csv ${outdir}/frac_r1 0 $start_allele_1 $end_allele_1 $start_allele_2 $end_allele_2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

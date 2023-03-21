process REPEAT_DIST_DISTANCE {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val allele_number
    val outdir

    output:
    path "*/stat_r1/*.stat.csv",   emit: stat_raw
    path "*/stat_r2/*.stat.csv",   emit: stat_r2
    path "*/plot_r1/*.png",        emit: plot_r1
    path "*/plot_r2/*.png",        emit: plot_r2
    tuple val(meta), path("*/count_r1/*.csv"),       emit: count_r1
    tuple val(meta), path("*/count_r2/*.csv"),       emit: count_r2
    path "*/frac_r1/*.csv",        emit: frac_r1
    path "*/frac_r2/*.csv",        emit: frac_r2
    path  "versions.yml",                                        emit: versions

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
    mkdir -p ${outdir}/stat_r1 ${outdir}/stat_r2 ${outdir}/plot_r1 ${outdir}/plot_r2 ${outdir}/count_r1 ${outdir}/count_r2

    repeat_dist_distance.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} ${outdir} $args

    # calculate fractions
    mkdir ${outdir}/frac_r1 ${outdir}/frac_r2

    if [ $allele_number -eq 1 ]; then
      calculate_frac.py $prefix ${outdir}/stat_r1/${prefix}.stat.csv ${outdir}/frac_r1 $length_cutoff_1_low $length_cutoff_1_high
      calculate_frac.py $prefix ${outdir}/stat_r2/${prefix}.stat.csv ${outdir}/frac_r2 $length_cutoff_1_low $length_cutoff_1_high
    else if [ $allele_number -eq 2 ]; then
      calculate_frac_2.py $prefix ${outdir}/stat_r1/${prefix}.stat.csv ${outdir}/frac_r1 $length_cutoff_1_low $length_cutoff_1_high $length_cutoff_2_low $length_cutoff_2_high
      calculate_frac_2.py $prefix ${outdir}/stat_r2/${prefix}.stat.csv ${outdir}/frac_r2 $length_cutoff_1_low $length_cutoff_1_high $length_cutoff_2_low $length_cutoff_2_high
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

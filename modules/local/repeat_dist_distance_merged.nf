process REPEAT_DIST_DISTANCE_MERGED {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val allele_number
    val outdir

    output:
    // path "*/stat/*.csv",                    emit: stat
    path "*/stat/*.stat.csv",               emit: stat_raw
    path "*/plot/*.png",                    emit: plot
    tuple val(meta), path("*/count/*.csv"), emit: count
    path "*/frac/*.csv",                    emit: frac
    tuple val(meta), path("*/frac/*.csv"),  emit: frac_meta
    path  "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def length_cutoff_1_low  = "${meta.length_cutoff_1_low}"
    def length_cutoff_1_high = "${meta.length_cutoff_1_high}"
    def length_cutoff_2_low  = "${meta.length_cutoff_2_low}" ?: 1 // simply place holder
    def length_cutoff_2_high = "${meta.length_cutoff_2_high}" ?: 2

    // def args_frac = task.ext.args_frac ?: ''
    def prefix = "${meta.id}"

    """
    mkdir -p ${outdir}/stat ${outdir}/plot ${outdir}/count ${outdir}/frac

    repeat_dist_distance_merged.py ${prefix}.fastq.gz ${prefix} ${outdir} $args

    # calculate fractions
    if [ $allele_number -eq 1 ]; then
      calculate_frac.py $prefix ${outdir}/stat/${prefix}.stat.csv ${outdir}/frac $length_cutoff_1_low $length_cutoff_1_high
    elif [ $allele_number -eq 2 ]; then
      calculate_frac_2.py $prefix ${outdir}/stat/${prefix}.stat.csv ${outdir}/frac $length_cutoff_1_low $length_cutoff_1_high $length_cutoff_2_low $length_cutoff_2_high
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

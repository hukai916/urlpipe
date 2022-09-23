process REPEAT_DIST_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    val outdir

    output:
    path "*/cutoff_1/*",       emit: cutoff_1
    path "*/cutoff_3/*",       emit: cutoff_3
    path "*/cutoff_10/*",      emit: cutoff_10
    path "*/cutoff_30/*",      emit: cutoff_30
    path "*/cutoff_100/*",     emit: cutoff_100
    path "*/cutoff_1/frac/*.csv", emit: frac_1
    path "*/cutoff_3/frac/*.csv", emit: frac_3
    path "*/cutoff_10/frac/*.csv", emit: frac_10
    path "*/cutoff_30/frac/*.csv", emit: frac_30
    path "*/cutoff_100/frac/*.csv", emit: frac_100
    path  "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # plot repeat length distribution with UMI cutoff: 1, 3, 10, 30, 100

    mkdir -p ${outdir}/cutoff_1 ${outdir}/cutoff_3 ${outdir}/cutoff_10 ${outdir}/cutoff_30 ${outdir}/cutoff_100

    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_1 1 $args
    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_3 3 $args
    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_10 10 $args
    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_30 30 $args
    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_100 100 $args

    # calculate fractions
    # tem="$args_frac"
    # suffix="\${tem// /_}"
    #mkdir ${outdir}/cutoff_1/frac_\$suffix ${outdir}/cutoff_3/frac_\$suffix ${outdir}/cutoff_10/frac_\$suffix ${outdir}/cutoff_30/frac_\$suffix ${outdir}/cutoff_100/frac_\$suffix
    mkdir ${outdir}/cutoff_1/frac ${outdir}/cutoff_3/frac ${outdir}/cutoff_10/frac ${outdir}/cutoff_30/frac ${outdir}/cutoff_100/frac

    calculate_frac.py $prefix ${outdir}/cutoff_1/stat_mode_${prefix}_cutoff_1.csv ${outdir}/cutoff_1/frac "$args_frac"
    calculate_frac.py $prefix ${outdir}/cutoff_3/stat_mode_${prefix}_cutoff_3.csv ${outdir}/cutoff_3/frac "$args_frac"
    calculate_frac.py $prefix ${outdir}/cutoff_10/stat_mode_${prefix}_cutoff_10.csv ${outdir}/cutoff_10/frac "$args_frac"
    calculate_frac.py $prefix ${outdir}/cutoff_30/stat_mode_${prefix}_cutoff_30.csv ${outdir}/cutoff_30/frac "$args_frac"
    calculate_frac.py $prefix ${outdir}/cutoff_100/stat_mode_${prefix}_cutoff_100.csv ${outdir}/cutoff_100/frac "$args_frac"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

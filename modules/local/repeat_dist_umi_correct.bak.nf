process REPEAT_DIST_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    path stat_raw
    val outdir

    output:
    path "*/cutoff_1/*",       emit: cutoff_1 // not if "*/cutoff_1", the resume is problematic
    path "*/cutoff_3/*",       emit: cutoff_3
    path "*/cutoff_10/*",      emit: cutoff_10
    path "*/cutoff_30/*",      emit: cutoff_30
    path "*/cutoff_100/*",     emit: cutoff_100
    path "*/frac_1/*.csv",     emit: frac_1
    path "*/frac_3/*.csv",     emit: frac_3
    path "*/frac_10/*.csv",    emit: frac_10
    path "*/frac_30/*.csv",    emit: frac_30
    path "*/frac_100/*.csv",   emit: frac_100
    path "*/plot_std/*",             emit: plot // STD and Violin plot of read length after UMI correction.
    path stat_raw,             emit: stat_raw
    path "*/cutoff_1/stat_mode*.csv",       emit: cutoff_1_mode_stat
    path "*/cutoff_3/stat_mode*.csv",       emit: cutoff_3_mode_stat
    path "*/cutoff_10/stat_mode*.csv",      emit: cutoff_10_mode_stat
    path "*/cutoff_30/stat_mode*.csv",      emit: cutoff_30_mode_stat
    path "*/cutoff_100/stat_mode*.csv",     emit: cutoff_100_mode_stat
    // path "5d_r1_repeat_dist_umi_correct/cutoff_100/frac/*.csv", emit: frac_100 // this won't work, may cause a weird issue, could be Nextflow's problem.
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

    mkdir -p ${outdir}/frac_1 ${outdir}/frac_3 ${outdir}/frac_10 ${outdir}/frac_30 ${outdir}/frac_100

    mkdir -p ${outdir}/plot_std

    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_1 1 $args
    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_3 3 $args
    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_10 10 $args
    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_30 30 $args
    repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_100 100 $args

    # calculate fractions
    # tem="$args_frac"
    # suffix="\${tem// /_}"
    #mkdir 5d_r1_repeat_dist_umi_correct/cutoff_1/frac_\$suffix 5d_r1_repeat_dist_umi_correct/cutoff_3/frac_\$suffix 5d_r1_repeat_dist_umi_correct/cutoff_10/frac_\$suffix 5d_r1_repeat_dist_umi_correct/cutoff_30/frac_\$suffix 5d_r1_repeat_dist_umi_correct/cutoff_100/frac_\$suffix

    calculate_frac.py $prefix ${outdir}/cutoff_1/stat_mode_${prefix}_cutoff_1.csv ${outdir}/frac_1 "$args_frac"
    calculate_frac.py $prefix ${outdir}/cutoff_3/stat_mode_${prefix}_cutoff_3.csv ${outdir}/frac_3 "$args_frac"
    calculate_frac.py $prefix ${outdir}/cutoff_10/stat_mode_${prefix}_cutoff_10.csv ${outdir}/frac_10 "$args_frac"
    calculate_frac.py $prefix ${outdir}/cutoff_30/stat_mode_${prefix}_cutoff_30.csv ${outdir}/frac_30 "$args_frac"
    calculate_frac.py $prefix ${outdir}/cutoff_100/stat_mode_${prefix}_cutoff_100.csv ${outdir}/frac_100 "$args_frac"

    # plot sd of read lengths vs UMI cutoffs; as well as read length violin plot:
    plot_sd_read_length_umi_cutoff.py $prefix ${outdir}/plot_std $stat_raw ${outdir}/cutoff_1/stat_mode_${prefix}_cutoff_1.csv ${outdir}/cutoff_3/stat_mode_${prefix}_cutoff_3.csv ${outdir}/cutoff_10/stat_mode_${prefix}_cutoff_10.csv ${outdir}/cutoff_30/stat_mode_${prefix}_cutoff_30.csv ${outdir}/cutoff_100/stat_mode_${prefix}_cutoff_100.csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

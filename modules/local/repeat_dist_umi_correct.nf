process REPEAT_DIST_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    val outdir

    output:
    // path "5d_r1_repeat_dist_umi_correct/cutoff_1/*",       emit: cutoff_1
    // path "5d_r1_repeat_dist_umi_correct/cutoff_3/*",       emit: cutoff_3
    // path "5d_r1_repeat_dist_umi_correct/cutoff_10/*",      emit: cutoff_10
    // path "5d_r1_repeat_dist_umi_correct/cutoff_30/*",      emit: cutoff_30
    // path "5d_r1_repeat_dist_umi_correct/cutoff_100/*",     emit: cutoff_100
    path "5d_r1_repeat_dist_umi_correct/frac_1/*.csv", emit: frac_1
    path "5d_r1_repeat_dist_umi_correct/frac_3/*.csv", emit: frac_3
    path "5d_r1_repeat_dist_umi_correct/frac_10/*.csv", emit: frac_10
    path "5d_r1_repeat_dist_umi_correct/frac_30/*.csv", emit: frac_30
    path "5d_r1_repeat_dist_umi_correct/frac_100/*.csv", emit: frac_100
    path  "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # plot repeat length distribution with UMI cutoff: 1, 3, 10, 30, 100

    mkdir -p 5d_r1_repeat_dist_umi_correct/cutoff_1 5d_r1_repeat_dist_umi_correct/cutoff_3 5d_r1_repeat_dist_umi_correct/cutoff_10 5d_r1_repeat_dist_umi_correct/cutoff_30 5d_r1_repeat_dist_umi_correct/cutoff_100

    mkdir -p 5d_r1_repeat_dist_umi_correct/frac_1 5d_r1_repeat_dist_umi_correct/frac_3 5d_r1_repeat_dist_umi_correct/frac_10 5d_r1_repeat_dist_umi_correct/frac_30 5d_r1_repeat_dist_umi_correct/frac_100


    repeat_dist_umi_correct.py $csv $prefix 5d_r1_repeat_dist_umi_correct/cutoff_1 1 $args
    repeat_dist_umi_correct.py $csv $prefix 5d_r1_repeat_dist_umi_correct/cutoff_3 3 $args
    repeat_dist_umi_correct.py $csv $prefix 5d_r1_repeat_dist_umi_correct/cutoff_10 10 $args
    repeat_dist_umi_correct.py $csv $prefix 5d_r1_repeat_dist_umi_correct/cutoff_30 30 $args
    repeat_dist_umi_correct.py $csv $prefix 5d_r1_repeat_dist_umi_correct/cutoff_100 100 $args

    # calculate fractions
    # tem="$args_frac"
    # suffix="\${tem// /_}"
    #mkdir 5d_r1_repeat_dist_umi_correct/cutoff_1/frac_\$suffix 5d_r1_repeat_dist_umi_correct/cutoff_3/frac_\$suffix 5d_r1_repeat_dist_umi_correct/cutoff_10/frac_\$suffix 5d_r1_repeat_dist_umi_correct/cutoff_30/frac_\$suffix 5d_r1_repeat_dist_umi_correct/cutoff_100/frac_\$suffix

    calculate_frac.py $prefix 5d_r1_repeat_dist_umi_correct/cutoff_1/stat_mode_${prefix}_cutoff_1.csv 5d_r1_repeat_dist_umi_correct/frac_1 "$args_frac"
    calculate_frac.py $prefix 5d_r1_repeat_dist_umi_correct/cutoff_3/stat_mode_${prefix}_cutoff_3.csv 5d_r1_repeat_dist_umi_correct/frac_3 "$args_frac"
    calculate_frac.py $prefix 5d_r1_repeat_dist_umi_correct/cutoff_10/stat_mode_${prefix}_cutoff_10.csv 5d_r1_repeat_dist_umi_correct/frac_10 "$args_frac"
    calculate_frac.py $prefix 5d_r1_repeat_dist_umi_correct/cutoff_30/stat_mode_${prefix}_cutoff_30.csv 5d_r1_repeat_dist_umi_correct/frac_30 "$args_frac"
    calculate_frac.py $prefix 5d_r1_repeat_dist_umi_correct/cutoff_100/stat_mode_${prefix}_cutoff_100.csv 5d_r1_repeat_dist_umi_correct/frac_100 "$args_frac"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

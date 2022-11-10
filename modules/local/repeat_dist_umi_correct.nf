process REPEAT_DIST_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    path stat_raw
    val umi_cutoffs
    val outdir

    output:
    path "*/cutoff_*/*/*",                        emit: cutoff // not if "*/cutoff_1", the resume is problematic
    tuple val(meta), path("*/frac_*/mode/*.csv"), emit: frac_meta_mode
    path "*/frac_*/mode/*.csv",                   emit: frac_mode
    path "*/cutoff_*/mode/stat_mode*.csv",        emit: cutoff_mode_stat

    tuple val(meta), path("*/frac_*/ld/*.csv"),   emit: frac_meta_ld
    path "*/frac_*/ld/*.csv",                     emit: frac_ld
    path "*/cutoff_*/ld/stat_ld*.csv",            emit: cutoff_ld_stat

    path "*/plot_read_length_std*/*",             emit: plot // STD and Violin plot of read length after UMI correction.
    tuple val(meta), path("${outdir}/input/*.csv"),           emit: stat_raw_meta
    path "${outdir}/input/*.csv",                             emit: stat_raw
    // path "5d_r1_repeat_dist_umi_correct/cutoff_100/frac/*.csv", emit: frac_100 // this won't work, may cause a weird issue, could be Nextflow's problem.
    path  "versions.yml",      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ## Dev notes:
    ## remove all space from string: https://stackoverflow.com/questions/13659318/how-to-remove-space-from-string

    mkdir -p ${outdir}/input
    cp $stat_raw ${outdir}/input/

    ## PLOT1: plot repeat length distribution with UMI cutoff
    # calculate fractions
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/cutoff_\$i ${outdir}/frac_\$i
      mkdir ${outdir}/cutoff_\$i/mode ${outdir}/cutoff_\$i/mean ${outdir}/cutoff_\$i/ld
      mkdir ${outdir}/frac_\$i/mode ${outdir}/frac_\$i/mean ${outdir}/frac_\$i/ld

      repeat_dist_umi_correct.py $csv $prefix ${outdir}/cutoff_\$i \$i $args

      calculate_frac.py $prefix ${outdir}/cutoff_\$i/mode/stat_mode_${prefix}_cutoff_\$i.csv ${outdir}/frac_\$i/mode \$i "$args_frac"
      calculate_frac.py $prefix ${outdir}/cutoff_\$i/mean/stat_mean_${prefix}_cutoff_\$i.csv ${outdir}/frac_\$i/mean \$i "$args_frac"
      calculate_frac.py $prefix ${outdir}/cutoff_\$i/ld/stat_ld_${prefix}_cutoff_\$i.csv ${outdir}/frac_\$i/ld \$i "$args_frac"

    done

    ## PLOT2: plot std of read lengths vs UMI cutoffs
    mkdir -p ${outdir}/plot_read_length_std_mode
    stat_mode_csv=""
    for i in "\${umi_cutoffs_array[@]}"
    do
      stat_mode_csv+=" ${outdir}/cutoff_\$i/mode/stat_mode_${prefix}_cutoff_\$i.csv"
    done

    plot_sd_read_length_umi_cutoff.py $prefix ${outdir}/plot_read_length_std_mode "$umi_cutoffs" $stat_raw \$stat_mode_csv

    # for ld:
    mkdir -p ${outdir}/plot_read_length_std_ld
    stat_ld_csv=""
    for i in "\${umi_cutoffs_array[@]}"
    do
      stat_ld_csv+=" ${outdir}/cutoff_\$i/ld/stat_ld_${prefix}_cutoff_\$i.csv"
    done

    plot_sd_read_length_umi_cutoff.py $prefix ${outdir}/plot_read_length_std_ld "$umi_cutoffs" $stat_raw \$stat_ld_csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

process REPEAT_DIST_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    path stat_raw
    val umi_cutoffs
    val allele_number
    val outdir

    output:
    path "*/frac_above_below/cutoff_*/ld/*.csv",                     emit: frac_ld
    path "*/frac_above_below/cutoff_*/mode/*.csv",                   emit: frac_mode
    tuple val(meta), path("*/frac_above_below/cutoff_*/ld/*.csv"),   emit: frac_meta_ld
    tuple val(meta), path("*/frac_above_below/cutoff_*/mode/*.csv"), emit: frac_meta_mode

    path "*/read_length_distribution/cutoff_*/*/*",                  emit: cutoff // if "*/cutoff_1", the resume is problematic
    path "*/read_length_distribution/cutoff_*/mode/stat_mode*.csv",  emit: cutoff_mode_stat
    path "*/read_length_distribution/cutoff_*/ld/stat_ld*.csv",      emit: cutoff_ld_stat

    path "*/plot_along_cutoffs/plot_read_length_std*/*",             emit: plot // STD and Violin plot of read length after UMI correction.
    tuple val(meta), path("${outdir}/input/*.csv"),                  emit: stat_raw_meta
    path "${outdir}/input/*.csv",                                    emit: stat_raw
    // path "5d_r1_repeat_dist_umi_correct/cutoff_100/frac/*.csv", emit: frac_100 // this won't work, may cause a weird issue, could be Nextflow's problem.
    path  "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def length_cutoff_1_low  = "${meta.length_cutoff_1_low}"
    def length_cutoff_1_high = "${meta.length_cutoff_1_high}"
    def length_cutoff_2_low  = "${meta.length_cutoff_2_low}" ?: 1 // simply place holder
    def length_cutoff_2_high = "${meta.length_cutoff_2_high}" ?: 2
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${outdir}/input
    cp $stat_raw ${outdir}/input/

    # Plot1: repeat length distribution with different UMI cutoffs
      # calculate fractions
    umi_cutoffs_str="0,$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/read_length_distribution/cutoff_\$i ${outdir}/frac_above_below/cutoff_\$i
      mkdir ${outdir}/read_length_distribution/cutoff_\$i/mode ${outdir}/read_length_distribution/cutoff_\$i/mean ${outdir}/read_length_distribution/cutoff_\$i/ld
      mkdir ${outdir}/frac_above_below/cutoff_\$i/mode ${outdir}/frac_above_below/cutoff_\$i/mean ${outdir}/frac_above_below/cutoff_\$i/ld

      repeat_dist_umi_correct.py $csv $prefix ${outdir}/read_length_distribution/cutoff_\$i \$i $args

      if [ $allele_number -eq 1 ]; then
        calculate_frac.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/mode/stat_mode_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/cutoff_\$i/mode $length_cutoff_1_low $length_cutoff_1_high
        calculate_frac.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/mean/stat_mean_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/cutoff_\$i/mean $length_cutoff_1_low $length_cutoff_1_high
        calculate_frac.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/ld/stat_ld_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/cutoff_\$i/ld $length_cutoff_1_low $length_cutoff_1_high
      elif [ $allele_number -eq 2 ]; then
        calculate_frac_2.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/mode/stat_mode_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/cutoff_\$i/mode $length_cutoff_1_low $length_cutoff_1_high $length_cutoff_2_low $length_cutoff_2_high
        calculate_frac_2.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/mean/stat_mean_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/cutoff_\$i/mean $length_cutoff_1_low $length_cutoff_1_high $length_cutoff_2_low $length_cutoff_2_high
        calculate_frac_2.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/ld/stat_ld_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/cutoff_\$i/ld $length_cutoff_1_low $length_cutoff_1_high $length_cutoff_2_low $length_cutoff_2_high
      fi
    done

    # Plot2: plot std of read lengths vs UMI cutoffs
    mkdir -p ${outdir}/plot_along_cutoffs/plot_read_length_std_mode
    stat_mode_csv=""
    for i in "\${umi_cutoffs_array[@]}"
    do
      stat_mode_csv+=" ${outdir}/read_length_distribution/cutoff_\$i/mode/stat_mode_${prefix}_cutoff_\$i.csv"
    done

    plot_sd_read_length_umi_cutoff.py $prefix ${outdir}/plot_along_cutoffs/plot_read_length_std_mode "$umi_cutoffs" $stat_raw \$stat_mode_csv

    # for ld:
    mkdir -p ${outdir}/plot_along_cutoffs/plot_read_length_std_ld
    stat_ld_csv=""
    for i in "\${umi_cutoffs_array[@]}"
    do
      stat_ld_csv+=" ${outdir}/read_length_distribution/cutoff_\$i/ld/stat_ld_${prefix}_cutoff_\$i.csv"
    done

    plot_sd_read_length_umi_cutoff.py $prefix ${outdir}/plot_along_cutoffs/plot_read_length_std_ld "$umi_cutoffs" $stat_raw \$stat_ld_csv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    ## Dev notes:
    ## remove all space from string: https://stackoverflow.com/questions/13659318/how-to-remove-space-from-string


    """
}

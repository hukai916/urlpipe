process INDEL_READS_UMI_CORRECT {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.4"

    input:
    tuple val(meta_5p), path(reads_5p)
    tuple val(meta_3p), path(reads_3p)
    val umi_cutoffs
    val outdir

    output:
    path "${outdir}/*.png", emit: plot
    path  "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta_5p}"

    """
    mkdir -p ${outdir}/5p_3p ${outdir}/cutoffs

    # merge 5p and 3p:
    zcat $reads_5p $reads_3p > ${outdir}/5p_3p/${prefix}.fastq.gz

    # get umi count table:


    # process per UMI cutoffs:
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/read_length_distribution/cutoff_\$i ${outdir}/frac_above_below/frac_above_below/frac_\$i
      mkdir ${outdir}/read_length_distribution/cutoff_\$i/mode ${outdir}/read_length_distribution/cutoff_\$i/mean ${outdir}/read_length_distribution/cutoff_\$i/ld
      mkdir ${outdir}/frac_above_below/frac_\$i/mode ${outdir}/frac_above_below/frac_above_below/frac_\$i/mean ${outdir}/frac_above_below/frac_above_below/frac_\$i/ld

      repeat_dist_umi_correct.py $csv $prefix ${outdir}/read_length_distribution/cutoff_\$i \$i $args

      calculate_frac.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/mode/stat_mode_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/frac_\$i/mode \$i "$args_frac"
      calculate_frac.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/mean/stat_mean_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/frac_\$i/mean \$i "$args_frac"
      calculate_frac.py $prefix ${outdir}/read_length_distribution/cutoff_\$i/ld/stat_ld_${prefix}_cutoff_\$i.csv ${outdir}/frac_above_below/frac_\$i/ld \$i "$args_frac"

    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

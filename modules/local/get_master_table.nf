process GET_MASTER_TABLE {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    // path frac_count_cutoff_0
    // path frac_count_cutoff_x
    // path indel_count_cutoff_x
    // val umi_cutoffs
    // val mode
    // val outdir

      path repeat_frac_csv // repeat_length_fraction_umi_x.csv
      path indel_csv // read_count_umi_cutoff_x.csv
      path stat_repeat_length_distribution // repeat_length_count_x_umi_x.csv
      val umi_cutoffs
      val allele_number
      val repeat_bins
      

    output:
      path "*.csv",         emit: csv
      path "versions.yml",  emit: versions

    when:
      task.ext.when == null || task.ext.when

    script:
      def args = task.ext.args ?: ''

      """
      umi_cutoffs_str="0,$umi_cutoffs"
      umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
      
      # master table 1: total counts, indel counts, allele counts/fraction
      for i in "\${umi_cutoffs_array[@]}"
      do
        get_master_table_allele.py repeat_length_fraction_umi_\${i}_*.csv read_count_umi_cutoff_\${i}.csv $allele_number master_table_allele_umi_\${i}.csv
      done

      # master table 2: total counts, indel counts, repeat bin counts/fraction 
      for i in "\${umi_cutoffs_array[@]}"
      do
        get_master_table_repeat_bin.py repeat_length_count_*_umi_\${i}.csv read_count_umi_cutoff_\${i}.csv "$repeat_bins" master_table_repeat_bin_umi_\${i}.csv
      done

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """
}

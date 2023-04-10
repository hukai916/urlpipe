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

      path repeat_frac_csv
      path indel_csv
      val umi_cutoffs
      val allele_number

    output:
      path "versions.yml",  emit: versions

    when:
      task.ext.when == null || task.ext.when

    script:
      def args = task.ext.args ?: ''

      """
      umi_cutoffs_str="0,$umi_cutoffs"
      umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
      for i in "\${umi_cutoffs_array[@]}"
      do
        get_master_table.py repeat_length_fraction_umi_\$i.csv read_count_umi_cutoff_\$i.csv $allele_number master_table_umi_\$i_${prefix}.csv
      done

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """
}

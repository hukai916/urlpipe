process GET_BIN_PLOT {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
      path csv
      val repeat_bins

    output:
      path "*.html",         emit: html
      path "versions.yml",   emit: versions

    when:
      task.ext.when == null || task.ext.when

    script:
      // def bins = task.ext.bins ?: ''
      // def use_ratio = task.ext.use_ratio ?: 'no'
      def use_repeat_unit_bp = task.ext.use_repeat_unit_bp ?: 'no'
      def repeat_unit_bp = task.ext.repeat_unit_bp ?: ''

      """
      # Use raw count:
      // for x in *.csv; do
      //   filename="\${x%.csv}".html
      //   get_bin_plot.py "$repeat_bins" $use_ratio $use_repeat_unit_bp $repeat_unit_bp \$x \$filename
      // done

      # Use ratio:


      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """
}

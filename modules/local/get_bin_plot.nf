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
      # Four plots will be created: 
        # using counts: without indel, with indel 
        # using ratio: without indel, with indel

      # Use raw count:
      for x in *.csv; do
        filename="\${x%.csv}".count.withoutIndel.html
        get_bin_plot.py "$repeat_bins" "no" $use_repeat_unit_bp $repeat_unit_bp \$x "no" \$filename
      done

      for x in *.csv; do
        filename="\${x%.csv}".count.withIndel.html
        get_bin_plot.py "$repeat_bins" "no" $use_repeat_unit_bp $repeat_unit_bp \$x "yes" \$filename
      done

      # Use count ratio:
      for x in *.csv; do
        filename="\${x%.csv}".ratio.withoutIndel.html
        get_bin_plot.py "$repeat_bins" "yes" $use_repeat_unit_bp $repeat_unit_bp \$x "no" \$filename
      done

      for x in *.csv; do
        filename="\${x%.csv}".ratio.withIndel.html
        get_bin_plot.py "$repeat_bins" "yes" $use_repeat_unit_bp $repeat_unit_bp \$x "yes" \$filename
      done

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """
}

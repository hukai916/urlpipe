process GET_BIN_PLOT {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
      path csv

    output:
      path "*.html",         emit: html
      path "versions.yml",   emit: versions

    when:
      task.ext.when == null || task.ext.when

    script:
      def args = task.ext.args ?: ''

      """
      touch test.html

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """
}

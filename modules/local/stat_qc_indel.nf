process STAT_QC_INDEL {
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(no_indel_fastq) // just be obtain meta
    path mean_qc
    path per_site_qc

    output:
    path "*.html",        emit: stat
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def per_site_range_to_display = task.ext.per_site_range_to_display ?: '1000:1200'

    """
    # Step1: obtain mean read quality distribution per read category
    plot_mean_qc.py $mean_qc plot_mean_qc_${prefix}.html

    # Step2: obtain per mapped site read quality distribution per read 
    plot_per_site_qc.py $per_site_qc $per_site_range_to_display

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS
    """
}

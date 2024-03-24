process PREP_REF {
    label 'process_low'
    container "hukai916/miniconda3_bio:0.4"

    input:

    output:
    path "*.fasta",       emit: ref
    path  "versions.yml", emit: versions

    script:
    def reference = file(task.ext.reference) ?: ''
    def repeat_start = task.ext.repeat_start ?: ''
    def repeat_end = task.ext.repeat_end ?: ''
    def repeat_unit = task.ext.repeat_unit ?: ''
    def repeat_range = task.ext.repeat_range ?: ''

    """
    prep_ref.py "$reference" $repeat_start $repeat_end $repeat_unit $repeat_range ref_repeat_range.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
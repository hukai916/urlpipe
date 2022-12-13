process READ_LENGTH_DIST {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    path "*/fastq",   emit: fastq
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_frac = task.ext.args_frac ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    read_length_dist.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${prefix} ${outdir} $args

    # add cutoff_0:
    mkdir -p ${outdir}/fastq/cutoff_0
    cp ${prefix}_1.fastq.gz ${outdir}/fastq/cutoff_0/
    cp ${prefix}_2.fastq.gz ${outdir}/fastq/cutoff_0/
    l=\$(zcat ${prefix}_1.fastq.gz | awk 'NR%4 == 0 {print}' | wc -l)
    echo ${prefix},\$l > ${outdir}/fastq/cutoff_0/${prefix}.csv

    # cat stat:
    umi_cutoffs_str="0,"
    umi_cutoffs_str+="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      cat ${outdir}/fastq/cutoff_\$i/*.csv > ${outdir}/fastq/cutoff_\$i/all_sample.csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

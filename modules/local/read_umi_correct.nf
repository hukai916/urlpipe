process READ_UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    path reads
    val umi_cutoffs
    val outdir

    output:
    path "*/stat_r2/*.stat.csv",   emit: stat_r2
    path "*/plot_r1/*.png",        emit: plot_r1
    path "*/plot_r2/*.png",        emit: plot_r2
    tuple val(meta), path("*/count_r1/*.csv"),       emit: count_r1
    tuple val(meta), path("*/count_r2/*.csv"),       emit: count_r2
    path  "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    read_umi_correct.py $csv ${prefix} ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${outdir}/fastq/ "$umi_cutoffs"

    # add cutoff_0:
    mkdir -p ${outdir}/fastq/cutoff_0
    cp ${prefix}_1.fastq.gz ${outdir}/fastq/cutoff_0/
    cp ${prefix}_2.fastq.gz ${outdir}/fastq/cutoff_0/
    l=\$(zcat ${prefix}_1.fastq.gz | awk 'NR%4 == 0 {print}' | wc -l)
    echo ${prefix},\$l > ${outdir}/fastq/cutoff_0/${prefix}.csv

    # cat stat:
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      cat ${outdir}/fastq/cutoff_\$i/ld/*.csv > ${outdir}/fastq/cutoff_\$i/ld/all_sample.csv
      cat ${outdir}/fastq/cutoff_\$i/mode/*.csv > ${outdir}/fastq/cutoff_\$i/mode/all_sample.csv
    done
    cat ${outdir}/fastq/cutoff_0/*.csv > ${outdir}/fastq/cutoff_0/all_sample.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

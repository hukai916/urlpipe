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

    path "*/fastq/**/*.fastq.gz",     emit: fastq
    path "*/count/*/ld/*.csv",        emit: count_ld
    path "*/count/*/mode/*.csv",      emit: count_mode
    path "versions.yml",         emit: versions

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

    # mv count.csv to count folder
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/count/cutoff_\$i/mode ${outdir}/count/cutoff_\$i/ld
      mv ${outdir}/fastq/cutoff_\$i/mode/*.csv ${outdir}/count/cutoff_\$i/mode/
      mv ${outdir}/fastq/cutoff_\$i/ld/*.csv ${outdir}/count/cutoff_\$i/ld/
    done
    # for easy cat_stat_umi:
    mkdir -p ${outdir}/count/cutoff_0/ld ${outdir}/count/cutoff_0/mode
    cp ${outdir}/fastq/cutoff_0/*.csv ${outdir}/count/cutoff_0/ld/${prefix}_cutoff_0.csv
    mv ${outdir}/fastq/cutoff_0/*.csv ${outdir}/count/cutoff_0/mode/${prefix}_cutoff_0.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

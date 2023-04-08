process READ_COUNT_PER_UMI_CUTOFF {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(reads)
    val umi_cutoffs

    output:
    tuple val(meta), path("reads_per_umi_*.csv"),  emit: csv
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    # 1. obtain master table: UMI,read_count
    reads_per_umi.py ${prefix}_1.fastq.gz reads_per_umi_${prefix}.csv

    # 2. obtain read count per different umi cutoff
    umi_cutoffs_str="0,$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str#[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      read_count_per_umi_cutoff.py reads_per_umi_${prefix}.csv \$i read_count_per_umi_${prefix}_umi_cutoff_\$i.csv
    done


    # obtain master table: UMI, read_count
    # read_count_umi_correct.py 

    #     repeat_length_per_read_default.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz raw_repeat_length_per_read_default_${prefix}.csv $args


    # read_umi_correct.py $csv ${prefix} ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz ${outdir}/fastq/ "$umi_cutoffs"

    # # add cutoff_0:
    # mkdir -p ${outdir}/fastq/cutoff_0
    # cp ${prefix}_1.fastq.gz ${outdir}/fastq/cutoff_0/
    # cp ${prefix}_2.fastq.gz ${outdir}/fastq/cutoff_0/
    # l=\$(zcat ${prefix}_1.fastq.gz | awk 'NR%4 == 0 {print}' | wc -l)
    # echo ${prefix},\$l > ${outdir}/fastq/cutoff_0/${prefix}.csv

    # # mv count.csv to count folder
    # umi_cutoffs_str="$umi_cutoffs"
    # umi_cutoffs_array=(\$(echo \${umi_cutoffs_str#[[:blank:]]/} | tr "," " "))
    # for i in "\${umi_cutoffs_array[@]}"
    # do
    #   mkdir -p ${outdir}/count/cutoff_\$i/mode ${outdir}/count/cutoff_\$i/ld
    #   mv ${outdir}/fastq/cutoff_\$i/mode/*.csv ${outdir}/count/cutoff_\$i/mode/
    #   mv ${outdir}/fastq/cutoff_\$i/ld/*.csv ${outdir}/count/cutoff_\$i/ld/
    # done
    # # for easy cat_stat_umi:
    # mkdir -p ${outdir}/count/cutoff_0/ld ${outdir}/count/cutoff_0/mode
    # cp ${outdir}/fastq/cutoff_0/*.csv ${outdir}/count/cutoff_0/ld/${prefix}_cutoff_0.csv
    # mv ${outdir}/fastq/cutoff_0/*.csv ${outdir}/count/cutoff_0/mode/${prefix}_cutoff_0.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python #g" )
    END_VERSIONS

    """
}

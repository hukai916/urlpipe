process READ_COUNT_PER_UMI_CUTOFF {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/bioinfo:0.1"

    input:
    tuple val(meta), path(reads)
    val umi_cutoffs

    output:
    tuple val(meta), path("read_count_*.csv"),  emit: csv
    path "read_count_*.csv",                    emit: csv_pure
    path "fastq/*.gz",                          emit: fastq
    path "versions.yml",                        emit: versions

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
      read_count_per_umi_cutoff.py reads_per_umi_${prefix}.csv \$i read_count_${prefix}_umi_cutoff_\$i.csv
    done

    # 3. obtain reads after umi correction
    mkdir fastq
    for i in "\${umi_cutoffs_array[@]}"
    do
      get_reads_umi_cutoff.py \$i reads_per_umi_${prefix}.csv ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz fastq/indel_${prefix}_umi_\${i}_R1.fastq fastq/indel_${prefix}_umi_\${i}_R2.fastq
    done
    gzip fastq/*.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python #g" )
    END_VERSIONS

    """
}

process INDEL_READS_UMI_CORRECT {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.4"

    input:
    tuple val(meta_5p), path(reads_5p)
    tuple val(meta_3p), path(reads_3p)
    val umi_cutoffs
    val outdir

    output:
    path "${outdir}/*.png", emit: plot
    path  "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta_5p}"

    """
    mkdir -p ${outdir}/5p_3p ${outdir}/read_length_distribution

    # merge 5p and 3p:
    zcat $reads_5p $reads_3p > ${outdir}/5p_3p/${prefix}.tem.fastq.gz
    scutls fastq -u -i ${outdir}/5p_3p/${prefix}.tem.fastq.gz -o ${outdir}/5p_3p/${prefix}.fastq.gz
    rm ${outdir}/5p_3p/${prefix}.tem.fastq.gz

    # get umi count tables:
      ## table1: seq_id,read_length
    read_length_dist.py ${outdir}/5p_3p/${prefix}.fastq.gz $prefix > ${outdir}/5p_3p/${prefix}.csv
      ## table2: read_length,freq
    cat ${outdir}/5p_3p/${prefix}.csv | cut -d "," -f 2 | sort | uniq -c | sort -k 2n | awk 'BEGIN {OFS=""} { print \$2","\$1}' > ${outdir}/5p_3p/${prefix}.stat.csv

    # get read counts and fastq:
      ## count_0:
    mkdir ${outdir}/read_length_distribution/cutoff_0
    cat ${outdir}/5p_3p/${prefix}.csv | wc -l > ${outdir}/read_length_distribution/cutoff_0/${prefix}.count.txt
    cp ${outdir}/5p_3p/${prefix}.fastq.gz ${outdir}/read_length_distribution/cutoff_0/

      ## count at different UMI cutoffs:
    umi_cutoffs_str="$umi_cutoffs"
    umi_cutoffs_array=(\$(echo \${umi_cutoffs_str//[[:blank:]]/} | tr "," " "))
    for i in "\${umi_cutoffs_array[@]}"
    do
      mkdir -p ${outdir}/cutoffs/cutoff_\$i
      length_dist_umi_correct.py ${outdir}/5p_3p/${prefix}.csv $prefix ${outdir}/read_length_distribution/cutoff_\$i \$i
    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}

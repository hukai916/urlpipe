// include non-repeat_stat_merge specific modules as XXX_MERGE in order not to be overwritten by other config files.
include { CLASSIFY_MERGE } from '../modules/local/classify_merge'
include { STAT_CSV_MERGE } from '../modules/local/stat_csv_merge'
include { FASTQC_SINGLE as FASTQC_SINGLE_MERGE } from '../modules/local/fastqc_single'
include { FASTQC as FASTQC_MERGE } from '../modules/nf-core/modules/fastqc/main'
include { REPEAT_LENGTH_DISTRIBUTION_MERGE } from '../modules/local/repeat_length_distribution_merge'
include { STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT as STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE } from '../modules/local/stat_repeat_length_distribution_default'
include { REPEAT_LENGTH_DISTRIBUTION_PER_UMI as REPEAT_LENGTH_DISTRIBUTION_PER_UMI_MERGE } from '../modules/local/repeat_length_distribution_per_umi'
include { PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI as PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI_MERGE } from '../modules/local/plot_repeat_length_distribution_per_umi'
include { REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT as REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT } from '../modules/local/repeat_length_distribution_default_umi_correct'
include { STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT as STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT } from '../modules/local/stat_repeat_length_distribution_default_umi_correct'
include { REPEAT_LENGTH_FRACTION as REPEAT_LENGTH_FRACTION_MERGE } from '../modules/local/repeat_length_fraction'

workflow REPEAT_STAT_MERGE {
    take:
      reads 
      ch_versions

    main:
      //
      // MODULE: BBmerge
      // 3_read_category/3d_classify_merge
      CLASSIFY_MERGE ( reads )
      ch_versions = ch_versions.mix(CLASSIFY_MERGE.out.versions)
      // 3_read_category/3d_classify_merge/classify_merge.csv
      STAT_CSV_MERGE (CLASSIFY_MERGE.out.csv.collect())

      //
      // MODULE: FastQC for merged and unmerged reads
      // 2_qc_and_umi/2a_fastqc/fastq_merge
      FASTQC_SINGLE_MERGE ( CLASSIFY_MERGE.out.reads )
      ch_versions = ch_versions.mix(FASTQC_SINGLE_MERGE.out.versions)
      // 2_qc_and_umi/2a_fastqc/fastq_non_merge
      FASTQC_MERGE ( CLASSIFY_MERGE.out.reads_others )
      ch_versions = ch_versions.mix(FASTQC_MERGE.out.versions)

      //
      // MODULE: repeat length distribution determined with merged reads
      // not published
      REPEAT_LENGTH_DISTRIBUTION_MERGE ( CLASSIFY_MERGE.out.reads )
      ch_versions = ch_versions.mix(REPEAT_LENGTH_DISTRIBUTION_MERGE.out.versions)
      // 4_repeat_statistics/4a_repeat_length_distribution/repeat_length_count_merge_umi_0.csv|html
      STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE (REPEAT_LENGTH_DISTRIBUTION_MERGE.out.repeat_length_count_merge_pure.collect())

      //
      // MODULE: repeat length count distribution per umi group
      // 4_repeat_statistics/4b_repeat_length_distribution_per_umi/csv
      REPEAT_LENGTH_DISTRIBUTION_PER_UMI_MERGE (
        REPEAT_LENGTH_DISTRIBUTION_MERGE.out.repeat_length_per_read_merge,
      params.umi_cutoffs )
      // 4_repeat_statistics/4b_repeat_length_distribution_per_umi/html
      PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI_MERGE (
        REPEAT_LENGTH_DISTRIBUTION_PER_UMI_MERGE.out.csv,
        params.umi_cutoffs )

      //
      // MODULE: repeat length count per umi group corrected
      // test_not_publish
      REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT (
        REPEAT_LENGTH_DISTRIBUTION_MERGE.out.repeat_length_per_read_merge,
        params.umi_correction_method,
        params.umi_cutoffs)
      // 4_repeat_statistics/4a_repeat_length_distribution/repeat_length_count_default_umi_x.csv|html
      STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT (REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT.out.repeat_length_count_default_umi_correct.collect(),
      params.umi_cutoffs)

      //
      //  MODULE: obtain fraction above and below for each sample at each cutoff
      // 4_repeat_statistics/4c_repeat_length_fraction
      REPEAT_LENGTH_FRACTION_MERGE (
        // just to obtain the sample meta info:
        REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT.out.umi_readcount_readlength_corrected,
        // master table for umi_0:
        STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE.out.csv,
        // master table for umi_x:
        STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT.out.csv,
        params.allele_number,
        params.umi_cutoffs,
        )

      stat_table = STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE.out.csv.mix(
        STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT.out.csv
      ).collect()

    emit:
      csv      = stat_table
      csv_frac = REPEAT_LENGTH_FRACTION_MERGE.out.csv.collect()
      versions = ch_versions
}

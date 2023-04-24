include { CLASSIFY_MERGE } from '../modules/local/classify_merge'
include { STAT_CSV_MERGE } from '../modules/local/stat_csv_merge'
include { FASTQC_SINGLE  } from '../modules/local/fastqc_single'
include { FASTQC         } from '../modules/nf-core/modules/fastqc/main'
include { REPEAT_LENGTH_DISTRIBUTION_MERGE } from '../modules/local/repeat_length_distribution_merge'
include { STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT as STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE } from '../modules/local/stat_repeat_length_distribution_default'
include { REPEAT_LENGTH_DISTRIBUTION_PER_UMI } from '../modules/local/repeat_length_distribution_per_umi'
include { PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI } from '../modules/local/plot_repeat_length_distribution_per_umi'
include { REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT as REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT } from '../modules/local/repeat_length_distribution_default_umi_correct'
include { STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT as STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT } from '../modules/local/stat_repeat_length_distribution_default_umi_correct'
include { REPEAT_LENGTH_FRACTION as REPEAT_LENGTH_FRACTION_TEST } from '../modules/local/repeat_length_fraction'

//
// if mode == "merge"
//

// still need to load FASTQC1 module even if it has been loaded in urlpipe.nf

// include { REPEAT_DIST_DISTANCE_MERGED       } from '../modules/local/repeat_dist_distance_merged'
// include { REPEAT_DIST_WITHIN_UMI_GROUP as REPEAT_DIST_WITHIN_UMI_GROUP_MERGE } from '../modules/local/repeat_dist_within_umi_group'
// include { UMI_GROUP_STAT as UMI_GROUP_STAT_MERGE } from '../modules/local/umi_group_stat'

// include { REPEAT_DIST_UMI_CORRECT as REPEAT_DIST_UMI_CORRECT_MERGE } from '../modules/local/repeat_dist_umi_correct'
// include { PLOT_FRAC_CUTOFF as PLOT_FRAC_CUTOFF_MERGE } from '../modules/local/plot_frac_cutoff'
// include { PLOT_FRAC_CUTOFF as PLOT_FRAC_CUTOFF_MERGE_2 } from '../modules/local/plot_frac_cutoff'
// include { PLOT_UMI_GROUPS as PLOT_UMI_GROUPS_MERGE   } from '../modules/local/plot_umi_groups'
// include { PLOT_UMI_GROUPS as PLOT_UMI_GROUPS_MERGE_2 } from '../modules/local/plot_umi_groups'
// include { PLOT_FRAC_UMI_CUTOFF as PLOT_FRAC_UMI_CUTOFF_MERGE } from '../modules/local/plot_frac_umi_cutoff'
// include { PLOT_FRAC_UMI_CUTOFF as PLOT_FRAC_UMI_CUTOFF_MERGE_2 } from '../modules/local/plot_frac_umi_cutoff'

// include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_MERGE }   from '../modules/local/cat_stat_cutoff'
// include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_MERGE_2}   from '../modules/local/cat_stat_cutoff'
// include { CAT_STAT as CAT_STAT5; CAT_STAT as CAT_STAT4 } from '../modules/local/cat_stat'

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
      FASTQC_SINGLE ( CLASSIFY_MERGE.out.reads )
      ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)
      // 2_qc_and_umi/2a_fastqc/fastq_non_merge
      FASTQC ( CLASSIFY_MERGE.out.reads_others )
      ch_versions = ch_versions.mix(FASTQC.out.versions)

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
      REPEAT_LENGTH_DISTRIBUTION_PER_UMI (
        REPEAT_LENGTH_DISTRIBUTION_MERGE.out.repeat_length_per_read_merge,
      params.umi_cutoffs )
      // 4_repeat_statistics/4b_repeat_length_distribution_per_umi/html
      PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI (
        REPEAT_LENGTH_DISTRIBUTION_PER_UMI.out.csv,
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
      // test_not_publish/4c
      REPEAT_LENGTH_FRACTION_TEST (
        // just to obtain the sample meta info:
        REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT.out.umi_readcount_readlength_corrected,
        // master table for umi_0:
        STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE.out.csv,
        // master table for umi_x:
        STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT.out.csv,
        params.allele_number,
        params.umi_cutoffs,
        )

      // stat_table = STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE.out.csv.mix(
      //   STAT_REPEAT_LENGTH_DISTRIBUTION_MERGE_UMI_CORRECT.out.csv
      // ).collect()




    // //
    // // MODULE: repeat distribution distance with merged reads
    // //
    // REPEAT_DIST_DISTANCE_MERGED (
    //   BBMERGE.out.reads_merged,
    //   params.allele_number,
    //   "4d_merge_repeat_distribution_distance"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_DISTANCE_MERGED.out.versions)

    // //
    // // MODULE: repeat distribution within UMI group
    // //
    // REPEAT_DIST_WITHIN_UMI_GROUP_MERGE (
    //   REPEAT_DIST_DISTANCE_MERGED.out.count,
    //   "5b_merge_repeat_dist_within_umi_group"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_WITHIN_UMI_GROUP_MERGE.out.versions)

    // //
    // // MODULE: UMI group statistics
    // //
    // UMI_GROUP_STAT_MERGE (
    //   REPEAT_DIST_DISTANCE_MERGED.out.count,
    //   REPEAT_DIST_DISTANCE_MERGED.out.stat_raw,
    //   "5c_merge_umi_group_stat"
    //   )
    // ch_versions = ch_versions.mix(UMI_GROUP_STAT_MERGE.out.versions)

    // //
    // // MODULE: repeat dist UMI corrected: 5d
    // //
    // REPEAT_DIST_UMI_CORRECT_MERGE (
    //   UMI_GROUP_STAT_MERGE.out.stat,
    //   UMI_GROUP_STAT_MERGE.out.stat_raw,
    //   params.umi_cutoffs,
    //   params.allele_number,
    //   "5d_merge_repeat_dist_umi_correct"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_UMI_CORRECT_MERGE.out.versions)

    // CAT_STAT5_MERGE (
    //   REPEAT_DIST_DISTANCE_MERGED.out.frac.collect(),
    //   "4d_merge_repeat_distribution_distance/frac_merge",
    //   "all_sample_frac",
    //   "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT5_MERGE.out.versions)

    // CAT_STAT_CUTOFF_MERGE ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_mode.collect(),
    //                         "mode",
    //                         "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std",
    //                         params.umi_cutoffs,
    //                         "all_sample_frac",
    //                         "5d_merge_repeat_dist_umi_correct/frac_above_below" ) // subdir: frac_xxx
    // CAT_STAT_CUTOFF_MERGE_2 ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_ld.collect(),
    // "ld", "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std",
    // params.umi_cutoffs,
    // "all_sample_frac",
    // "5d_merge_repeat_dist_umi_correct/frac_above_below" ) // subdir: frac_xxx
    // PLOT_FRAC_CUTOFF_MERGE (
    //   // CAT_STAT5_MERGE.out.stat,
    //   CAT_STAT_CUTOFF_MERGE.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw.collect(),
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_mode_stat.collect(),
    //   params.umi_cutoffs,
    //   "all_sample_frac",
    //   "5d_merge_repeat_dist_umi_correct/plot_all_samples/mode" // plot_read_length_violin and plot_frac_barplot
    // )
    // PLOT_FRAC_CUTOFF_MERGE_2 (
    //   // CAT_STAT5_MERGE.out.stat,
    //   CAT_STAT_CUTOFF_MERGE_2.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw.collect(),
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_ld_stat.collect(),
    //   params.umi_cutoffs,
    //   "all_sample_frac",
    //   "5d_merge_repeat_dist_umi_correct/plot_all_samples/ld" // plot_read_length_violin and plot_frac_barplot
    // )
    // PLOT_UMI_GROUPS_MERGE (
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw_meta,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_mode_stat,
    //   params.umi_cutoffs,
    //   "mode",
    //   "5d_merge_repeat_dist_umi_correct/plot_along_cutoffs/plot_umi_groups_mode"
    // )
    // PLOT_UMI_GROUPS_MERGE_2 (
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw_meta,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_ld_stat,
    //   params.umi_cutoffs,
    //   "ld",
    //   "5d_merge_repeat_dist_umi_correct/plot_along_cutoffs/plot_umi_groups_ld"
    // )
    // PLOT_FRAC_UMI_CUTOFF_MERGE (
    //   REPEAT_DIST_DISTANCE_MERGED.out.frac.collect(),
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_meta_mode,
    //   params.umi_cutoffs,
    //   "5d_merge_repeat_dist_umi_correct/plot_along_cutoffs/plot_frac_umi_cutoff_mode"
    // )
    // PLOT_FRAC_UMI_CUTOFF_MERGE_2 (
    //   REPEAT_DIST_DISTANCE_MERGED.out.frac.collect(),
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_meta_ld,
    //   params.umi_cutoffs,
    //   "5d_merge_repeat_dist_umi_correct/plot_along_cutoffs/plot_frac_umi_cutoff_ld"
    // )

    emit:
    versions = ch_versions
}

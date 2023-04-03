//
// if use_read == "merge"
//

// still need to load FASTQC1 module even if it has been loaded in urlpipe.nf
include { FASTQC_SINGLE               } from '../modules/local/fastqc_single'

include { BBMERGE                           } from '../modules/local/bbmerge'
include { REPEAT_DIST_DISTANCE_MERGED       } from '../modules/local/repeat_dist_distance_merged'
include { REPEAT_DIST_WITHIN_UMI_GROUP as REPEAT_DIST_WITHIN_UMI_GROUP_MERGE } from '../modules/local/repeat_dist_within_umi_group'
include { UMI_GROUP_STAT as UMI_GROUP_STAT_MERGE } from '../modules/local/umi_group_stat'

include { REPEAT_DIST_UMI_CORRECT as REPEAT_DIST_UMI_CORRECT_MERGE } from '../modules/local/repeat_dist_umi_correct'
include { PLOT_FRAC_CUTOFF as PLOT_FRAC_CUTOFF_MERGE } from '../modules/local/plot_frac_cutoff'
include { PLOT_FRAC_CUTOFF as PLOT_FRAC_CUTOFF_MERGE_2 } from '../modules/local/plot_frac_cutoff'
include { PLOT_UMI_GROUPS as PLOT_UMI_GROUPS_MERGE   } from '../modules/local/plot_umi_groups'
include { PLOT_UMI_GROUPS as PLOT_UMI_GROUPS_MERGE_2 } from '../modules/local/plot_umi_groups'
include { PLOT_FRAC_UMI_CUTOFF as PLOT_FRAC_UMI_CUTOFF_MERGE } from '../modules/local/plot_frac_umi_cutoff'
include { PLOT_FRAC_UMI_CUTOFF as PLOT_FRAC_UMI_CUTOFF_MERGE_2 } from '../modules/local/plot_frac_umi_cutoff'

include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_MERGE }   from '../modules/local/cat_stat_cutoff'
include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_MERGE_2}   from '../modules/local/cat_stat_cutoff'
include { CAT_STAT as CAT_STAT5; CAT_STAT as CAT_STAT4 } from '../modules/local/cat_stat'

workflow USE_READ_MERGE {
    take:
    reads_through
    ch_versions

    main:
    //
    // MODULE: BBmerge
    //
    BBMERGE ( reads_through )
    ch_versions = ch_versions.mix(BBMERGE.out.versions)

    //
    // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    //
    CAT_STAT4 (
      BBMERGE.out.stat.collect(),
      "4b_bbmerge/stat",
      "all_sample",
      "sample_name,count_merge,count_non_merge" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT4.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC_SINGLE (
      BBMERGE.out.reads_merged
      // "4c_merge_fastqc"
      )
    ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)

    //
    // MODULE: repeat distribution distance with merged reads
    //
    REPEAT_DIST_DISTANCE_MERGED (
      BBMERGE.out.reads_merged,
      params.allele_number,
      "4d_merge_repeat_distribution_distance"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_DISTANCE_MERGED.out.versions)

    //
    // MODULE: repeat distribution within UMI group
    //
    REPEAT_DIST_WITHIN_UMI_GROUP_MERGE (
      REPEAT_DIST_DISTANCE_MERGED.out.count,
      "5b_merge_repeat_dist_within_umi_group"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_WITHIN_UMI_GROUP_MERGE.out.versions)

    //
    // MODULE: UMI group statistics
    //
    UMI_GROUP_STAT_MERGE (
      REPEAT_DIST_DISTANCE_MERGED.out.count,
      REPEAT_DIST_DISTANCE_MERGED.out.stat_raw,
      "5c_merge_umi_group_stat"
      )
    ch_versions = ch_versions.mix(UMI_GROUP_STAT_MERGE.out.versions)

    //
    // MODULE: repeat dist UMI corrected: 5d
    //
    REPEAT_DIST_UMI_CORRECT_MERGE (
      UMI_GROUP_STAT_MERGE.out.stat,
      UMI_GROUP_STAT_MERGE.out.stat_raw,
      params.umi_cutoffs,
      params.allele_number,
      "5d_merge_repeat_dist_umi_correct"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_UMI_CORRECT_MERGE.out.versions)

    CAT_STAT5_MERGE (
      REPEAT_DIST_DISTANCE_MERGED.out.frac.collect(),
      "4d_merge_repeat_distribution_distance/frac_merge",
      "all_sample_frac",
      "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT5_MERGE.out.versions)

    CAT_STAT_CUTOFF_MERGE ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_mode.collect(),
                            "mode",
                            "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std",
                            params.umi_cutoffs,
                            "all_sample_frac",
                            "5d_merge_repeat_dist_umi_correct/frac_above_below" ) // subdir: frac_xxx
    CAT_STAT_CUTOFF_MERGE_2 ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_ld.collect(),
    "ld", "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std",
    params.umi_cutoffs,
    "all_sample_frac",
    "5d_merge_repeat_dist_umi_correct/frac_above_below" ) // subdir: frac_xxx
    PLOT_FRAC_CUTOFF_MERGE (
      // CAT_STAT5_MERGE.out.stat,
      CAT_STAT_CUTOFF_MERGE.out.stat,
      REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw.collect(),
      REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_mode_stat.collect(),
      params.umi_cutoffs,
      "all_sample_frac",
      "5d_merge_repeat_dist_umi_correct/plot_all_samples/mode" // plot_read_length_violin and plot_frac_barplot
    )
    PLOT_FRAC_CUTOFF_MERGE_2 (
      // CAT_STAT5_MERGE.out.stat,
      CAT_STAT_CUTOFF_MERGE_2.out.stat,
      REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw.collect(),
      REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_ld_stat.collect(),
      params.umi_cutoffs,
      "all_sample_frac",
      "5d_merge_repeat_dist_umi_correct/plot_all_samples/ld" // plot_read_length_violin and plot_frac_barplot
    )
    PLOT_UMI_GROUPS_MERGE (
      REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw_meta,
      REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_mode_stat,
      params.umi_cutoffs,
      "mode",
      "5d_merge_repeat_dist_umi_correct/plot_along_cutoffs/plot_umi_groups_mode"
    )
    PLOT_UMI_GROUPS_MERGE_2 (
      REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw_meta,
      REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_ld_stat,
      params.umi_cutoffs,
      "ld",
      "5d_merge_repeat_dist_umi_correct/plot_along_cutoffs/plot_umi_groups_ld"
    )
    PLOT_FRAC_UMI_CUTOFF_MERGE (
      REPEAT_DIST_DISTANCE_MERGED.out.frac.collect(),
      REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_meta_mode,
      params.umi_cutoffs,
      "5d_merge_repeat_dist_umi_correct/plot_along_cutoffs/plot_frac_umi_cutoff_mode"
    )
    PLOT_FRAC_UMI_CUTOFF_MERGE_2 (
      REPEAT_DIST_DISTANCE_MERGED.out.frac.collect(),
      REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_meta_ld,
      params.umi_cutoffs,
      "5d_merge_repeat_dist_umi_correct/plot_along_cutoffs/plot_frac_umi_cutoff_ld"
    )

    emit:
    versions = ch_versions
}

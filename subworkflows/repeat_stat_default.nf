//
// if params,use_read == "R1"
//

// include { REPEAT_LENGTH_DISTRIBUTION_DEFAULT     } from '../modules/local/REPEAT_LENGTH_DISTRIBUTION_DEFAULT'
include { REPEAT_LENGTH_DISTRIBUTION_DEFAULT  } from '../modules/local/repeat_length_distribution_default'
include { STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT    } from '../modules/local/stat_repeat_length_distribution_default'
include { REPEAT_LENGTH_DISTRIBUTION_PER_UMI  } from '../modules/local/repeat_length_distribution_per_umi'
include { PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI  } from '../modules/local/plot_repeat_length_distribution_per_umi'
include { REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT } from '../modules/local/repeat_length_distribution_default_umi_correct'



include { REPEAT_DIST_UMI_CORRECT as REPEAT_DIST_UMI_CORRECT_R1 } from '../modules/local/repeat_dist_umi_correct'

include { REPEAT_DIST_WITHIN_UMI_GROUP as REPEAT_DIST_WITHIN_UMI_GROUP_R1 } from '../modules/local/repeat_dist_within_umi_group'

include { UMI_GROUP_STAT } from '../modules/local/umi_group_stat'


// include { PLOT_FRAC as PLOT_FRAC_4D_R1    } from '../modules/local/plot_frac'
// include { PLOT_FRAC as PLOT_FRAC_4D_MERGE } from '../modules/local/plot_frac'
include { PLOT_FRAC_CUTOFF as PLOT_FRAC_CUTOFF_R1    } from '../modules/local/plot_frac_cutoff'

include { PLOT_FRAC_CUTOFF as PLOT_FRAC_CUTOFF_R1_2    } from '../modules/local/plot_frac_cutoff'

include { PLOT_UMI_GROUPS as PLOT_UMI_GROUPS_R1      } from '../modules/local/plot_umi_groups'
include { PLOT_UMI_GROUPS as PLOT_UMI_GROUPS_R1_2    } from '../modules/local/plot_umi_groups'

include { PLOT_FRAC_UMI_CUTOFF as PLOT_FRAC_UMI_CUTOFF_R1    } from '../modules/local/plot_frac_umi_cutoff'

include { PLOT_FRAC_UMI_CUTOFF as PLOT_FRAC_UMI_CUTOFF_R1_2  } from '../modules/local/plot_frac_umi_cutoff'

include { CAT_STAT as CAT_STAT5 } from '../modules/local/cat_stat'
include { CAT_STAT_CUTOFF             }   from '../modules/local/cat_stat_cutoff'
include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_2      }   from '../modules/local/cat_stat_cutoff'

workflow REPEAT_STAT_DEFAULT {
    take:
      reads
      ch_versions

    main:
    //
    // MODULE: repeat length distribution determined with R1/R2 reads
    // not published
    REPEAT_LENGTH_DISTRIBUTION_DEFAULT ( reads )
    ch_versions = ch_versions.mix(REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.versions)
    // 4_repeat_statistics/repeat_length_count_default_umi_0.csv|html
    STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT (REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.repeat_length_count_default_pure.collect())

    //
    // MODULE: repeat length count distribution per umi group
    // 4_repeat_statistics/4b_repeat_length_distribution_per_umi/csv
    REPEAT_LENGTH_DISTRIBUTION_PER_UMI (
      REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.repeat_length_per_read_default,
    params.umi_cutoffs )
    // 4_repeat_statistics/4b_repeat_length_distribution_per_umi/html
    PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI (
      REPEAT_LENGTH_DISTRIBUTION_PER_UMI.out.csv,
      params.umi_cutoffs )

    //
    // MODULE: repeat length count per umi group corrected
    // test_not_publish
    REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT (
      REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.repeat_length_per_read_default,
      params.umi_correction_method,
      params.umi_cutoffs)
    // 4_repeat_statistics/repeat_length_count_default_umi_x.csv|html
    // STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT (REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT.out.repeat_length_count_default_umi_correct.collect(),
    // params.umi_cutoffs)


// 4a_repeat_length_distribution

    //
    // MODULE: repeat length per umi group
    // 4_repeat_statistics/4b_repeat_length_per_umi/




    // // 2
    // // MODULE: UMI group stat: UMI read_count mean mode: 5c
    // //
    // UMI_GROUP_STAT (
    //   REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.count_r1,
    //   REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.stat_raw, // stat_raw store the raw stat before UMI correction
    //   "5c_r1_umi_group_stat"
    //   )
    // ch_versions = ch_versions.mix(UMI_GROUP_STAT.out.versions)

    // // 3
    // // MODULE: repeat dist UMI corrected: 5d
    // //
    // REPEAT_DIST_UMI_CORRECT_R1 (
    //   UMI_GROUP_STAT.out.stat,
    //   UMI_GROUP_STAT.out.stat_raw,
    //   params.umi_cutoffs,
    //   params.allele_number,
    //   "5d_r1_repeat_dist_umi_correct"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_UMI_CORRECT_R1.out.versions)

    // // 4
    // // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    // //
    // // REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.frac_r1.collect().view()
    // CAT_STAT5 (
    //   REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.frac_r1.collect(),
    //   "4d_r1_repeat_distribution_distance/frac_r1",
    //   "all_sample_frac",
    //   "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT5.out.versions)

    //
    // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    //
    // CAT_STAT6 (
    //   REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.frac_r2.collect(),
    //   "4d_repeat_distribution_distance/frac_r2",
    //   "all_sample_frac",
    //   "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT6.out.versions)


    // // // 5
    // // // MODULE: combine REPEAT_DIST_UMI_CORRECT_R1.out.frac_x into one file
    // // //
    // // // REPEAT_DIST_UMI_CORRECT_R1.out.frac_1.collect().view()
    // // CAT_STAT_CUTOFF is specifically designed for mode
    // CAT_STAT_CUTOFF ( REPEAT_DIST_UMI_CORRECT_R1.out.frac_mode.collect(),
    //                   "mode",
    //                   "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std",
    //                   params.umi_cutoffs,
    //                   "all_sample_frac",
    //                   "5d_r1_repeat_dist_umi_correct/frac_above_below" ) // subdir: frac_xxx

    // // For ld 6
    // CAT_STAT_CUTOFF_2 ( REPEAT_DIST_UMI_CORRECT_R1.out.frac_ld.collect(),
    // "ld", "sample_name,below_count,below_frac,below_mean,below_std,between_count,between_frac,beetween_mean,beetween_std,above_count,above_frac,above_mean,above_std",
    // params.umi_cutoffs,
    // "all_sample_frac",
    // "5d_r1_repeat_dist_umi_correct/frac_above_below" ) // subdir: frac_xxx

    // MODULE: plot the mean, std, and frac of all_sample.csv for frac stat.
    // Below has be integrated into PLOT_FRAC_CUTOFF_R1
    // PLOT_FRAC_4D_R1 (
    //   CAT_STAT5.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_R1.out.stat_raw.collect(),
    //   Channel.value(0),
    //   "4d_repeat_distribution_distance/plot_r1_frac"
    // )
    //
    // PLOT_FRAC_4D_MERGE (
    //   CAT_STAT5_MERGE.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw.collect(),
    //   Channel.value(0),
    //   "4d_merge_repeat_distribution_distance/plot_merge_frac"
    // )

    // add cutoff_0: without UMI correction, while
    // for mode:
    // PLOT_FRAC_CUTOFF_R1 (
    //   // CAT_STAT5.out.stat, // all_sample_frac_cutoff_0.csv
    //   CAT_STAT_CUTOFF.out.stat, // all_sample_frac_cutoff_x.csv
    //   REPEAT_DIST_UMI_CORRECT_R1.out.stat_raw.collect(), // individual stat without UMI correction
    //   REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_mode_stat.collect(),
    //   params.umi_cutoffs,
    //   "all_sample_frac",
    //   "5d_r1_repeat_dist_umi_correct/plot_all_samples/mode" // plot_read_length_violin and plot_frac_barplot
    // )
    //
    // // for ld:
    // PLOT_FRAC_CUTOFF_R1_2 (
    //   // CAT_STAT5.out.stat,
    //   CAT_STAT_CUTOFF_2.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_R1.out.stat_raw.collect(),
    //   REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_ld_stat.collect(),
    //   params.umi_cutoffs,
    //   "all_sample_frac",
    //   "5d_r1_repeat_dist_umi_correct/plot_all_samples/ld" // plot_read_length_violin and plot_frac_barplot
    // )
    //
    // // MODULE: plot UMI groups at different UMI cutoffs
    // PLOT_UMI_GROUPS_R1 (
    //   REPEAT_DIST_UMI_CORRECT_R1.out.stat_raw_meta,
    //   REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_mode_stat,
    //   params.umi_cutoffs,
    //   "mode",
    //   "5d_r1_repeat_dist_umi_correct/plot_along_cutoffs/plot_umi_groups_mode"
    // )
    //
    // // for ld:
    // PLOT_UMI_GROUPS_R1_2 (
    //   REPEAT_DIST_UMI_CORRECT_R1.out.stat_raw_meta,
    //   REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_ld_stat,
    //   params.umi_cutoffs,
    //   "ld",
    //   "5d_r1_repeat_dist_umi_correct/plot_along_cutoffs/plot_umi_groups_ld"
    // )
    //
    // // MODULE: plot above/below fraction at different UMI cutoffs
    // PLOT_FRAC_UMI_CUTOFF_R1 (
    //   REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.frac_r1.collect(),
    //   REPEAT_DIST_UMI_CORRECT_R1.out.frac_meta_mode,
    //   params.umi_cutoffs,
    //   "5d_r1_repeat_dist_umi_correct/plot_along_cutoffs/plot_frac_umi_cutoff_mode"
    // )
    //
    // // for ld:
    // PLOT_FRAC_UMI_CUTOFF_R1_2 (
    //   REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.frac_r1.collect(),
    //   REPEAT_DIST_UMI_CORRECT_R1.out.frac_meta_ld,
    //   params.umi_cutoffs,
    //   "5d_r1_repeat_dist_umi_correct/plot_along_cutoffs/plot_frac_umi_cutoff_ld"
    // )

    emit:
    // stat5         = CAT_STAT5.out.stat
    // cutoff_stat   = CAT_STAT_CUTOFF.out.stat
    // cutoff2_stat  = CAT_STAT_CUTOFF_2.out.stat
    versions      = ch_versions
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

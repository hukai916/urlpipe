include { READ_UMI_CORRECT } from '../modules/local/read_umi_correct'
include { READ_LENGTH_DIST } from '../modules/local/read_length_dist'
include { REPEAT_DIST_UMI_CORRECT as REPEAT_DIST_UMI_CORRECT_INDEL } from '../modules/local/repeat_dist_umi_correct'
include { UMI_GROUP_STAT as UMI_GROUP_STAT_INDEL } from '../modules/local/umi_group_stat'
include { COUNT_SUMMARY as COUNT_SUMMARY_LD                  } from '../modules/local/count_summary'
include { COUNT_SUMMARY as COUNT_SUMMARY_MODE                } from '../modules/local/count_summary'
include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_INDEL   } from '../modules/local/cat_stat_cutoff'
include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_INDEL_2 } from '../modules/local/cat_stat_cutoff'

include { READ_COUNT_PER_UMI_CUTOFF                  } from '../modules/local/read_count_per_umi_cutoff'
include { STAT_READ_COUNT_PER_UMI_CUTOFF             } from '../modules/local/stat_read_count_per_umi_cutoff'

workflow INDEL_STAT {
    take:
      reads
      reads_pure
      ch_versions

    main:

    // 
    // MODULE: read count umi correct
    // 5_indel_stat/5a_read_count_per_umi_cutoff/raw_csv
    READ_COUNT_PER_UMI_CUTOFF (
      reads,
      params.umi_cutoffs
    )

    // 
    // MODULE: stat read count umi correct
    // 
    STAT_READ_COUNT_PER_UMI_CUTOFF (
      READ_COUNT_PER_UMI_CUTOFF.out.csv.collect(),
      params.umi_cutoffs
    )

    // MODULE: INDEL reads distribution:
    READ_LENGTH_DIST (
      reads,
      "XXX_4d_indel_read_length_distribution"
    )
    ch_versions = ch_versions.mix(READ_LENGTH_DIST.out.versions)

    // MODULE: INDEL reads UMI stat
    UMI_GROUP_STAT_INDEL (
      READ_LENGTH_DIST.out.count_r1,
      READ_LENGTH_DIST.out.stat_raw, // stat_raw store the raw stat before UMI correction
      "XXX_5c_indel_umi_group_stat"
      )
    ch_versions = ch_versions.mix(UMI_GROUP_STAT_INDEL.out.versions)

    // MODULE: UMI correct
    // // if UMI cutoff filter results to 0 counts, the following module hicups, need more debugging
    // REPEAT_DIST_UMI_CORRECT_INDEL (
    //   UMI_GROUP_STAT_INDEL.out.stat,
    //   UMI_GROUP_STAT_INDEL.out.stat_raw,
    //   params.umi_cutoffs,
    //   "5d_indel_read_length_dist_umi_correct"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_UMI_CORRECT_INDEL.out.versions)

    // MODULE: INDEL UMI correct
    READ_UMI_CORRECT (
      UMI_GROUP_STAT_INDEL.out.stat,
      reads_pure.collect(),
      params.umi_cutoffs,
      "XXX_5e_indel_read_umi_correct"
      )
    ch_versions = ch_versions.mix(READ_UMI_CORRECT.out.versions)

    CAT_STAT_CUTOFF_INDEL ( READ_UMI_CORRECT.out.count_ld.collect(), "ld", "sample_name,read_count",
    "0," + params.umi_cutoffs,
    "all_sample_indel",
    "XXX_5e_indel_read_umi_correct/count" )

    CAT_STAT_CUTOFF_INDEL_2 ( READ_UMI_CORRECT.out.count_mode.collect(), "mode", "sample_name,read_count",
    "0," + params.umi_cutoffs,
    "all_sample_indel",
    "XXX_5e_indel_read_umi_correct/count" )

    emit:
    versions      = ch_versions

}

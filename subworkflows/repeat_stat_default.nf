include { REPEAT_LENGTH_DISTRIBUTION_DEFAULT  } from '../modules/local/repeat_length_distribution_default'
include { PREP_REF } from '../modules/local/prep_ref'
include { BWA } from '../modules/local/bwa'
include { BWA_LENGTH } from '../modules/local/bwa_length'
include { STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT    } from '../modules/local/stat_repeat_length_distribution_default'
include { REPEAT_LENGTH_DISTRIBUTION_PER_UMI  } from '../modules/local/repeat_length_distribution_per_umi'
include { PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI  } from '../modules/local/plot_repeat_length_distribution_per_umi'
include { REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT } from '../modules/local/repeat_length_distribution_default_umi_correct'
include { STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT    } from '../modules/local/stat_repeat_length_distribution_default_umi_correct'
include { REPEAT_LENGTH_FRACTION } from '../modules/local/repeat_length_fraction'

workflow REPEAT_STAT_DEFAULT {
    take:
      reads
      length_mode // how to determine the length of the repeat: distance-count based or reference-alignment based
      ch_versions

    main:
      //
      // MODULE: repeat length distribution determined with R1/R2 reads
      // not published

      if (length_mode == "distance_count") {
        REPEAT_LENGTH_DISTRIBUTION_DEFAULT ( reads )

        REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.repeat_length_per_read_default.set( {repeat_length_per_read_default} )
        REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.repeat_length_count_default_pure.set( {repeat_length_count_default_pure} )
        // ch_versions = ch_versions.mix(REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.versions)
      } else if (length_mode == "reference_align") {
        PREP_REF () // prepare references of varying repeat length
        BWA (PREP_REF.out.ref, reads)
        BWA_LENGTH (BWA.out.bam_reads)

        BWA_LENGTH.out.repeat_length_per_read_default.set( {repeat_length_per_read_default} )
        BWA_LENGTH.out.repeat_length_count_default_pure.set( {repeat_length_count_default_pure} )
      }

      // 4_repeat_statistics/4a_repeat_length_distribution/repeat_length_count_default_umi_0.csv|html
      STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT (repeat_length_count_default_pure.collect())

      //
      // MODULE: repeat length count distribution per umi group
      // 4_repeat_statistics/4b_repeat_length_distribution_per_umi/csv
      REPEAT_LENGTH_DISTRIBUTION_PER_UMI ( repeat_length_per_read_default, params.umi_cutoffs )

      // 4_repeat_statistics/4b_repeat_length_distribution_per_umi/html
      PLOT_REPEAT_LENGTH_DISTRIBUTION_PER_UMI (
        REPEAT_LENGTH_DISTRIBUTION_PER_UMI.out.csv,
        params.umi_cutoffs )

      //
      // MODULE: repeat length count per umi group corrected
      // test_not_publish
      REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT (
        repeat_length_per_read_default,
        params.umi_correction_method,
        params.umi_cutoffs)

      // 4_repeat_statistics/4a_repeat_length_distribution/repeat_length_count_default_umi_x.csv|html
      STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT (
        REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT.out.repeat_length_count_default_umi_correct.collect(),
        params.umi_cutoffs)

      //
      //  MODULE: obtain fraction above and below for each sample at each cutoff
      // test_not_publish/4c
      REPEAT_LENGTH_FRACTION (
        // just to obtain the sample meta info:
        REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT.out.umi_readcount_readlength_corrected,
        // master table for umi_0:
        STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.csv,
        // master table for umi_x:
        STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT.out.csv,
        params.allele_number,
        params.umi_cutoffs,
        )

      stat_table = STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.csv.mix(
        STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT.out.csv
      ).collect()

    emit:
      csv      = stat_table
      csv_frac = REPEAT_LENGTH_FRACTION.out.csv.collect()
      stat_repeat_length_distribution = STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT_UMI_CORRECT.out.csv.mix(
        STAT_REPEAT_LENGTH_DISTRIBUTION_DEFAULT.out.csv
      ).collect()
      versions = ch_versions
}

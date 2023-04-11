include { READ_COUNT_PER_UMI_CUTOFF      } from '../modules/local/read_count_per_umi_cutoff'
include { STAT_READ_COUNT_PER_UMI_CUTOFF } from '../modules/local/stat_read_count_per_umi_cutoff'

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
      // 5_indel_stat/5a_read_count_per_umi_cutoff/
      STAT_READ_COUNT_PER_UMI_CUTOFF (
        READ_COUNT_PER_UMI_CUTOFF.out.csv_pure.collect(),
        params.umi_cutoffs
      )

    emit:
      csv      = STAT_READ_COUNT_PER_UMI_CUTOFF.out.csv 
      versions = ch_versions
}

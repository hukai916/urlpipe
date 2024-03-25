include { GET_MASTER_TABLE } from '../modules/local/get_master_table'
include { GET_BIN_PLOT } from '../modules/local/get_bin_plot'

workflow GET_SUMMARY {
    take:
      repeat_frac_csv
      indel_csv
      umi_cutoffs
      allele_number
      stat_repeat_length_distribution
      ch_version

    main:
      //
      // MODULE: get_master_table
      // 6_summary/6a_master_table
      GET_MASTER_TABLE (repeat_frac_csv, indel_csv, umi_cutoffs, allele_number)

      // MODULE: get_bin_plot
      // 6_summary/6b_bin_plot
      GET_BIN_PLOT (stat_repeat_length_distribution)

    emit:
      versions = ch_version  
}

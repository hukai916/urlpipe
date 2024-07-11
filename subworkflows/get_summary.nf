include { GET_MASTER_TABLE } from '../modules/local/get_master_table'
include { GET_BIN_PLOT } from '../modules/local/get_bin_plot'

workflow GET_SUMMARY {
    take:
      repeat_frac_csv // used for alleles, calculated separately since we allow different start/end for alleles for different samples
      indel_csv
      umi_cutoffs
      allele_number // used for binning purpose
      repeat_bins // used for binning purpose with more flexibility
      stat_repeat_length_distribution // repeat_length, sample_1, sample2, xxx
      ch_version

    main:
      //
      // MODULE: get_master_table
      // 6_summary/6a_master_table
      GET_MASTER_TABLE (repeat_frac_csv, indel_csv, stat_repeat_length_distribution, umi_cutoffs, allele_number, repeat_bins)

      // MODULE: get_bin_plot
      // 6_summary/6b_bin_plot
      // stat_repeat_length_distribution.view()
      GET_BIN_PLOT (stat_repeat_length_distribution, repeat_bins)

    emit:
      versions = ch_version  
}

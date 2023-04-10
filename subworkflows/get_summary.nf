include { GET_MASTER_TABLE } from '../modules/local/get_master_table'


workflow GET_SUMMARY {
    take:
      repeat_csv
      indel_csv
      umi_cutoffs
      allele_number
      ch_version

    main:
      //
      // MODULE: get_master_table
      // 6_summary/6a_master_table
      GET_MASTER_TABLE (repeat_csv, indel_csv, umi_cutoffs, allele_number)

    emit:
      versions = ch_versions  
}

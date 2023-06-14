include { ADAPTOR_COUNT as ADAPTOR_COUNT_AP1    } from '../modules/local/adaptor_count'
include { ADAPTOR_COUNT as ADAPTOR_COUNT_AP2    } from '../modules/local/adaptor_count'
include { ADAPTOR_COUNT as ADAPTOR_COUNT_AP3_1  } from '../modules/local/adaptor_count'
include { ADAPTOR_COUNT as ADAPTOR_COUNT_AP3_2  } from '../modules/local/adaptor_count'
include { ADAPTOR_COUNT as ADAPTOR_COUNT_AP4_1  } from '../modules/local/adaptor_count'
include { ADAPTOR_COUNT as ADAPTOR_COUNT_AP4_2  } from '../modules/local/adaptor_count'
include { ADAPTOR_COUNT as ADAPTOR_COUNT_AP5    } from '../modules/local/adaptor_count'
// include { STAT_ADAPTOR } from '../modules/local/stat_ADAPTOR'

workflow ADAPTOR_COUNT_WF {
    take:
      reads

    main:
      // 
      // MODULE: ADAPTOR_COUNT
      // 1_preprocess_nanopore/1b_adaptor_count
      ADAPTOR_COUNT_AP1 ( reads, "raw_reads_AP1" )
      ADAPTOR_COUNT_AP2 ( reads, "raw_reads_AP2" )
      ADAPTOR_COUNT_AP3_1 ( reads, "raw_reads_AP3_1" )
      ADAPTOR_COUNT_AP3_2 ( reads, "raw_reads_AP3_2" )
      ADAPTOR_COUNT_AP4_1 ( reads, "raw_reads_AP4_1" )
      ADAPTOR_COUNT_AP4_2 ( reads, "raw_reads_AP4_2" )
      ADAPTOR_COUNT_AP5 ( reads, "raw_reads_AP5" )

      // 
      // MODULE: merge all stat into a single csv file 
      // 1_preprocess_nanopore/1b_adaptor_count
      ch_csv = ADAPTOR_COUNT_AP1.out.count_csv
                .mix( ADAPTOR_COUNT_AP2.out.count_csv )
                .mix( ADAPTOR_COUNT_AP3_1.out.count_csv )
                .mix( ADAPTOR_COUNT_AP3_2.out.count_csv )
                .mix( ADAPTOR_COUNT_AP4_1.out.count_csv )
                .mix( ADAPTOR_COUNT_AP4_2.out.count_csv )
                .mix( ADAPTOR_COUNT_AP5.out.count_csv )
                .collect()

    emit:
      csv      = ch_csv
      versions = ADAPTOR_COUNT_AP1.out.versions
}

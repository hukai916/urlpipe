include { BARCODE_COUNT as BARCODE_COUNT_BC1    } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC2    } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC3_1  } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC3_2  } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC4_1  } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC4_2  } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC5    } from '../modules/local/barcode_count'
// include { STAT_BARCODE } from '../modules/local/stat_barcode'

workflow BARCODE_COUNT_WF {
    take:
      reads

    main:
      // 
      // MODULE: BARCODE_COUNT
      // 1_preprocess_nanopore/1b_barcode_count
      BARCODE_COUNT_BC1 ( reads, "raw_reads_bc1" )
      BARCODE_COUNT_BC2 ( reads, "raw_reads_bc2" )
      BARCODE_COUNT_BC3_1 ( reads, "raw_reads_bc3_1" )
      BARCODE_COUNT_BC3_2 ( reads, "raw_reads_bc3_2" )
      BARCODE_COUNT_BC4_1 ( reads, "raw_reads_bc4_1" )
      BARCODE_COUNT_BC4_2 ( reads, "raw_reads_bc4_2" )
      BARCODE_COUNT_BC5 ( reads, "raw_reads_bc5" )

      // 
      // MODULE: merge all stat into a single csv file 
      // 1_preprocess_nanopore/1b_barcode_count
      ch_csv = BARCODE_COUNT_BC1.out.count_csv
                .mix( BARCODE_COUNT_BC2.out.count_csv )
                .mix( BARCODE_COUNT_BC3_1.out.count_csv )
                .mix( BARCODE_COUNT_BC3_2.out.count_csv )
                .mix( BARCODE_COUNT_BC4_1.out.count_csv )
                .mix( BARCODE_COUNT_BC4_2.out.count_csv )
                .mix( BARCODE_COUNT_BC5.out.count_csv )
                .collect()
      // STAT_BARCODE( ch_csv )

    emit:
      csv      = ch_csv,
      versions = BARCODE_COUNT_BC1.out.versions
}

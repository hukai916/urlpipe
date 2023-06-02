include { BARCODE_COUNT as BARCODE_COUNT_BC1    } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC2    } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC3_1  } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC3_2  } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC4_1  } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC4_2  } from '../modules/local/barcode_count'
include { BARCODE_COUNT as BARCODE_COUNT_BC5    } from '../modules/local/barcode_count'


workflow BARCODE_COUNT_WF {
    take:
      reads

    main:
      // 
      // MODULE: BARCODE_COUNT
      // 1_preprocess_nanopore/1b_barcode_count_bc1
      BARCODE_COUNT_BC1 ( reads, "raw_reads_bc1" )
      BARCODE_COUNT_BC2 ( reads, "raw_reads_bc2" )

    // emit:
    //   // reads    = ch_out_reads
    //   versions = ch_versions
}

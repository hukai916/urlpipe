include { BARCODE_COUNT_WF                     } from '../subworkflows/barcode_count_wf'
include { STAT_BARCODE                         } from '../modules/local/stat_barcode'
include { GET_VALID_NANOPORE_READS             } from '../modules/local/get_valid_nanopore_reads'
include { STAT_BARCODE as STAT_VALID_READS     } from '../modules/local/stat_barcode'

// to process forward reads:
include { CUTADAPT as CUTADAPT_NANOPORE_BC01   } from '../modules/nf-core/modules/cutadapt/main'
// include { DEMULPLEX }                            from '../modules/local/demultiplex'


include { CUTADAPT as CUTADAPT_NANOPORE_3END   } from '../modules/nf-core/modules/cutadapt/main'

include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { UMI_EXTRACT                 } from '../modules/local/umi_extract'
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/modules/fastqc/main'
include { FASTQC as FASTQC_CUTADAPT   } from '../modules/nf-core/modules/fastqc/main'
include { READ_PER_UMI                } from '../modules/local/read_per_umi'

workflow PREPROCESS_NANOPORE {
    take:
      reads
      ch_versions

    main:
      // 
      // MODULE: BARCODE_COUNT: using raw reads
      // 1_preprocess_nanopore/1a_barcode_count
      BARCODE_COUNT_WF ( reads )
      ch_versions = ch_versions.mix(BARCODE_COUNT_WF.out.versions)
      STAT_BARCODE ( reads, BARCODE_COUNT_WF.out.csv )

      // 
      // MODULE: GET_VALID_READS: valid means that read must contain both bc1 and bc2 in a row or rc
      // 1_preprocess_nanopore/1b_valid_reads
      GET_VALID_NANOPORE_READS ( reads )
      STAT_VALID_READS ( reads, GET_VALID_NANOPORE_READS.out.csv )

      // 
      // MODULE: CUTADAPT: bc01
      // 1_preprocess_nanopore/1c_cutadapt_bc01
      CUTADAPT_NANOPORE_BC01 ( GET_VALID_NANOPORE_READS.out.reads_valid )

      // 
      // MODULE: DEMULTIPLEX: using bc02
      // 1_preprocess_nanopore/1d_demultiplex
      // DEMULPLEX ( CUTADAPT_NANOPORE_BC01.out.reads )

      // For reads with forward direction:
      // TRIM 5END (BC1)
      // DEMULPLEX using BC2
      // TRIM 5END again (BC3)
      // UMI_EXTRACT
      // TRIM 5END (BC4) and 4END (BC5)
      // GET_SOLID_READS: MODULE: both 5' 200bp and 3' 200bp are mapped to hs or mm
      // got to DOWNSTREAM

      // For reads with reverse direction: 5' and 3' all relative as forward reads
      // TRIM 5END (reads_rc's 3END)
      // DEMULTIPLEX using BC2

      // 
      // MODULE: CUTADAPT_NANOPORE_5END
      // 1_preprocess_nanopore/1b_cutadapt_5end
      // CUTADAPT_NANOPORE_5END ( reads )
      // ch_versions = ch_versions.mix(CUTADAPT_NANOPORE_5END.out.versions)
     
      // 
      // MODULE: CUTADAPT_NANOPORE_3END
      // 1_preprocess_nanopore/1b_cutadapt_3end
      // CUTADAPT_NANOPORE_3END ( CUTADAPT_NANOPORE_5END.out.reads )




      //
      // MODULE: BARCODE_COUNT_TRIM
      // 1_preprocess_nanopore/1b_barcode_count_trim


      //
      // MODULE: FASTQC_RAW
      // 2_qc_and_umi/2a_fastqc/01_raw
      // FASTQC_RAW ( reads )
      // ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)

      // //
      // // MODULE: Cat Fastq
      // // 1_preprocess/1a_lane_merge
      // CAT_FASTQ ( reads )
      // ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

      // if (mode != "nanopore") {
      //   //
      //   // MODULE: UMI extract
      //   // 1_preprocess/1b_umi_extract
      //   UMI_EXTRACT ( CAT_FASTQ.out.reads )
      //   ch_versions = ch_versions.mix(UMI_EXTRACT.out.versions)

      //   //
      //   // MODULE: cutadapt
      //   // 1_preprocess/1c_cutadapt
      //   CUTADAPT ( UMI_EXTRACT.out.reads )
      //   ch_versions = ch_versions.mix(CUTADAPT.out.versions)
      //   ch_out_reads = CUTADAPT.out.reads

      //   //
      //   // MODULE: FastQC
      //   // 2_qc_and_umi/2a_fastqc/02_after_cutadapt
      //   FASTQC_CUTADAPT ( CUTADAPT.out.reads )
      //   ch_versions = ch_versions.mix(FASTQC_CUTADAPT.out.versions)

      //   //
      //   // MODULE: READ_PER_UMI: use after cutadapt reads
      //   // 2_qc_and_umi/read_per_umi_cutadapt
      //   READ_PER_UMI ( CUTADAPT.out.reads )
      //   ch_versions = ch_versions.mix(READ_PER_UMI.out.versions)
      // } else {
      //   ch_out_reads = CAT_FASTQ.out.reads
      // }


    emit:
      // reads    = ch_out_reads
      versions = ch_versions
}

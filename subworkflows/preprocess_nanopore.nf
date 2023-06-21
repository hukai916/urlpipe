include { ADAPTOR_COUNT_WF                     } from '../subworkflows/adaptor_count_wf'
include { STAT_ADAPTOR                         } from '../modules/local/stat_adaptor'
include { GET_VALID_NANOPORE_READS             } from '../modules/local/get_valid_nanopore_reads'
include { STAT_ADAPTOR as STAT_VALID_READS     } from '../modules/local/stat_adaptor'
include { CUTADAPT as CUTADAPT_NANOPORE_AP01   } from '../modules/nf-core/modules/cutadapt/main'
include { DEMULTIPLEX                          } from '../modules/local/demultiplex'
include { CUTADAPT_FASTQS as CUTADAPT_FASTQS_AP02 } from '../modules/nf-core/modules/cutadapt_fastqs/main'
include { UMI_EXTRACT_FASTQS                   } from '../modules/local/umi_extract_fastq'                   
include { CUTADAPT_FASTQS as CUTADAPT_FASTQS_AP03 } from '../modules/nf-core/modules/cutadapt_fastqs/main'
include { CUTADAPT_FASTQS as CUTADAPT_FASTQS_AP04 } from '../modules/nf-core/modules/cutadapt_fastqs/main'
include { GET_FULL_LENGTH_READS                } from '../modules/local/get_full_length_reads'
include { STAT as STAT_FULL_LENGTH             } from '../modules/local/stat'
include { SPLIT_ALLELE                         } from '../modules/local/split_allele'
include { STAT as STAT_SPLIT_ALLELE            } from '../modules/local/stat'




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
      // MODULE: ADAPTOR_COUNT: using raw reads
      // 1_preprocess_nanopore/1a_adaptor_count
      ADAPTOR_COUNT_WF ( reads )
      ch_versions = ch_versions.mix(ADAPTOR_COUNT_WF.out.versions)
      STAT_ADAPTOR ( reads, ADAPTOR_COUNT_WF.out.csv )

      // 
      // MODULE: GET_VALID_READS: valid means that read must contain both AP1 and AP2 in a row or rc
      // 1_preprocess_nanopore/1b_valid_reads
      GET_VALID_NANOPORE_READS ( reads )
      STAT_VALID_READS ( reads, GET_VALID_NANOPORE_READS.out.csv )

// 
// FOR forward reads combined with rc of rc_reads:
// 
      // 
      // MODULE: CUTADAPT: AP01
      // 1_preprocess_nanopore/1c_cutadapt_AP01
      CUTADAPT_NANOPORE_AP01 ( GET_VALID_NANOPORE_READS.out.reads_valid_combine )
      // CUTADAPT_NANOPORE_AP01_RC ( GET_VALID_NANOPORE_READS.out.reads_valid_rc )

      // 
      // MODULE: DEMULTIPLEX: using sample index (not AP02)
      // 1_preprocess_nanopore/1d_demultiplex/reads
      DEMULTIPLEX ( CUTADAPT_NANOPORE_AP01.out.reads, "0", "0::24" )

      // 
      // MODULE: CUTADAPT: AP02: AP02 is right after sample index
      // 1_preprocess_nanopore/1e_cutadapt_AP02
      CUTADAPT_FASTQS_AP02 ( DEMULTIPLEX.out.reads.flatten() )
      // flatten to spit split output fastq file individually to leverage pipeline parallelism

      // 
      // MODULE: UMI_EXTRACT_FASTQ
      // 1_preprocess_nanopore/1f_umi_extract
      UMI_EXTRACT_FASTQS ( CUTADAPT_FASTQS_AP02.out.reads )

      // 
      // MODULE: CUTADAPT: AP03: 2 versions; and AP04: 2 versions
      // 1_preprocess_nanopore/1g_cutadapt_ap03
      // 1_preprocess_nanopore/1g_cutadapt_ap04
      CUTADAPT_FASTQS_AP03 ( UMI_EXTRACT_FASTQS.out.reads )
      CUTADAPT_FASTQS_AP04 ( CUTADAPT_FASTQS_AP03.out.reads )
      
      // 
      // MODULE: GET_FULL_LENGTH_READS: start of ref fasta must appear in the beginning of the read and end of ref fasta must appear in the end of the read
      // 1_preprocess_nanopore/1h_full_length_read
      GET_FULL_LENGTH_READS ( CUTADAPT_FASTQS_AP04.out.reads, file(params.ref) )
      STAT_FULL_LENGTH ( GET_FULL_LENGTH_READS.out.stat.collect() )

      // 
      // MODULE: SPIT_ALLELE: split sample according to allele SNP information
      // 1_preprocess_nanopore/1i_split_allele
      if (params.allele_number == 2) {
        SPLIT_ALLELE ( GET_FULL_LENGTH_READS.out.reads )
        STAT_SPLIT_ALLELE ( SPLIT_ALLELE.out.stat.collect() )
      } 




      // For reads with forward direction:
      // TRIM 5END (AP1)
      // DEMULPLEX using AP2
      // TRIM 5END again (AP3)
      // UMI_EXTRACT
      // TRIM 5END (AP4) and 4END (AP5)
      // GET_SOLID_READS: MODULE: both 5' 200bp and 3' 200bp are mapped to hs or mm
      // got to DOWNSTREAM

      // For reads with reverse direction: 5' and 3' all relative as forward reads
      // TRIM 5END (reads_rc's 3END)
      // DEMULTIPLEX using AP2

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
      // MODULE: adaptor_COUNT_TRIM
      // 1_preprocess_nanopore/1b_adaptor_count_trim


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

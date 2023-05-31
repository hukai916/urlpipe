include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { UMI_EXTRACT                 } from '../modules/local/umi_extract'
include { CUTADAPT                    } from '../modules/nf-core/modules/cutadapt/main'
include { FASTQC as FASTQC_RAW        } from '../modules/nf-core/modules/fastqc/main'
include { FASTQC as FASTQC_CUTADAPT   } from '../modules/nf-core/modules/fastqc/main'
include { READ_PER_UMI                } from '../modules/local/read_per_umi'


workflow PREPROCESS_QC {
    take:
      reads
      mode
      ch_versions

    main:
      //
      // MODULE: FASTQC_RAW
      // 2_qc_and_umi/2a_fastqc/01_raw
      FASTQC_RAW ( reads )
      ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)

      //
      // MODULE: Cat Fastq
      // 1_preprocess/1a_lane_merge
      CAT_FASTQ ( reads )
      ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

      if (mode != "nanopore") {
        //
        // MODULE: UMI extract
        // 1_preprocess/1b_umi_extract
        UMI_EXTRACT ( CAT_FASTQ.out.reads )
        ch_versions = ch_versions.mix(UMI_EXTRACT.out.versions)

        //
        // MODULE: cutadapt
        // 1_preprocess/1c_cutadapt
        CUTADAPT ( UMI_EXTRACT.out.reads )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)
        ch_out_reads = CUTADAPT.out.reads

        //
        // MODULE: FastQC
        // 2_qc_and_umi/2a_fastqc/02_after_cutadapt
        FASTQC_CUTADAPT ( CUTADAPT.out.reads )
        ch_versions = ch_versions.mix(FASTQC_CUTADAPT.out.versions)

        //
        // MODULE: READ_PER_UMI: use after cutadapt reads
        // 2_qc_and_umi/read_per_umi_cutadapt
        READ_PER_UMI ( CUTADAPT.out.reads )
        ch_versions = ch_versions.mix(READ_PER_UMI.out.versions)
      } else {
        ch_out_reads = CAT_FASTQ.out.reads
      }


    emit:
      reads    = ch_out_reads
      versions = ch_versions
}

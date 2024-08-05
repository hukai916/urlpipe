include { CLASSIFY_LOCUS              } from '../modules/local/classify_locus'
include { STAT as STAT_LOCUS          } from '../modules/local/stat'
include { CLASSIFY_INDEL              } from '../modules/local/classify_indel'
include { STAT as STAT_INDEL          } from '../modules/local/stat'
include { CLASSIFY_READTHROUGH        } from '../modules/local/classify_readthrough'
include { STAT as STAT_READTHROUGH    } from '../modules/local/stat'
include { FASTQC as FASTQC_READTHROUGH} from '../modules/nf-core/modules/fastqc/main'
include { READ_PER_UMI as READ_PER_UMI_READTHROUGH } from '../modules/local/read_per_umi'

workflow CLASSIFY_READ {
    take:
      reads
      mode
      ch_versions

    main:
      //
      // MODULE: classify locus and stat
      // 3_read_category/3a_classify_locus
      CLASSIFY_LOCUS ( reads, file(params.ref) )
      ch_versions = ch_versions.mix(CLASSIFY_LOCUS.out.versions)
      STAT_LOCUS ( CLASSIFY_LOCUS.out.stat.collect() )
      ch_versions = ch_versions.mix(STAT_LOCUS.out.versions)

      //
      // MODULE: classify INDEL and stat
      // 3_read_category/3b_classify_indel
      CLASSIFY_INDEL ( CLASSIFY_LOCUS.out.reads_locus, file(params.ref))
      ch_versions = ch_versions.mix(CLASSIFY_INDEL.out.versions)
      STAT_INDEL ( CLASSIFY_INDEL.out.stat.collect() )
      ch_versions = ch_versions.mix(STAT_INDEL.out.versions)

      //
      // MODULE: classify_readthrough
      // 3_read_category/3c_classify_readthrough
      CLASSIFY_READTHROUGH ( CLASSIFY_INDEL.out.reads_no_indel )
      ch_versions = ch_versions.mix(CLASSIFY_READTHROUGH.out.versions)
      STAT_READTHROUGH ( CLASSIFY_READTHROUGH.out.stat.collect() )
      ch_versions = ch_versions.mix(STAT_READTHROUGH.out.versions)
      FASTQC_READTHROUGH ( reads )
      ch_versions = ch_versions.mix(FASTQC_READTHROUGH.out.versions)
      // 2_qc_and_umi/2c_read_per_umi_readthrough
      READ_PER_UMI_READTHROUGH ( CLASSIFY_READTHROUGH.out.reads_through )
      ch_versions = ch_versions.mix(READ_PER_UMI_READTHROUGH.out.versions)

      // MODULE: classify reads into alleles based on SNP
      // TODO module
      // if params.allele_number == 2:
      //   CLASSIFY_ALLELE ()

    emit:
      reads_indel_5p_or_3p      = CLASSIFY_INDEL.out.reads_indel_5p_or_3p
      reads_indel_5p_or_3p_pure = CLASSIFY_INDEL.out.reads_indel_5p_or_3p_pure
      reads_through = CLASSIFY_READTHROUGH.out.reads_through
      versions      = ch_versions
}

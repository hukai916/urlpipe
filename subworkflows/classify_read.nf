include { CLASSIFY_LOCUS              } from '../modules/local/classify_locus'
include { STAT as STAT_LOCUS          } from '../modules/local/stat'
include { CLASSIFY_INDEL              } from '../modules/local/classify_indel'
include { STAT as STAT_INDEL          } from '../modules/local/stat'
include { CLASSIFY_READTHROUGH        } from '../modules/local/classify_readthrough'
include { STAT as STAT_READTHROUGH    } from '../modules/local/stat'


workflow CLASSIFY_READ {
    take:
      reads
      ch_versions

    main:
      //
      // MODULE: classify locus and stat
      // 3_read_category/3a_classify_locus
      CLASSIFY_LOCUS ( reads )
      ch_versions = ch_versions.mix(CLASSIFY_LOCUS.out.versions)
      STAT_LOCUS ( CLASSIFY_LOCUS.out.stat.collect() )
      ch_versions = ch_versions.mix(STAT_LOCUS.out.versions)

      //
      // MODULE: classify INDEL and stat
      // 3_read_category/3b_classify_indel
      CLASSIFY_INDEL ( CLASSIFY_LOCUS.out.reads_locus )
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

    emit:
      reads_indel_5p_or_3p      = CLASSIFY_INDEL.out.reads_indel_5p_or_3p
      reads_indel_5p_or_3p_pure = CLASSIFY_INDEL.out.reads_indel_5p_or_3p_pure
      reads_through = CLASSIFY_READTHROUGH.out.reads_through
      ch_versions   = ch_versions
}

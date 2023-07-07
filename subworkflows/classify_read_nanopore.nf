include { CLASSIFY_LOCUS              } from '../modules/local/classify_locus'
include { STAT as STAT_LOCUS          } from '../modules/local/stat'
include { CLASSIFY_INDEL_NANOPORE     } from '../modules/local/classify_indel_nanopore'
include { STAT as STAT_INDEL          } from '../modules/local/stat'
include { STAT_QC_INDEL               } from '../modules/local/stat_qc_indel'

include { CLASSIFY_READTHROUGH        } from '../modules/local/classify_readthrough'
include { STAT as STAT_READTHROUGH    } from '../modules/local/stat'
include { FASTQC as FASTQC_READTHROUGH} from '../modules/nf-core/modules/fastqc/main'
include { READ_PER_UMI as READ_PER_UMI_READTHROUGH } from '../modules/local/read_per_umi'

workflow CLASSIFY_READ_NANOPORE {
    take:
      reads
      mode
      ch_versions

    main:
      //
      // MODULE: classify locus and stat
      // 3_read_category/3a_classify_locus
      // CLASSIFY_LOCUS ( reads )
      // ch_versions = ch_versions.mix(CLASSIFY_LOCUS.out.versions)
      // STAT_LOCUS ( CLASSIFY_LOCUS.out.stat.collect() )
      // ch_versions = ch_versions.mix(STAT_LOCUS.out.versions)

      //
      // MODULE: classify INDEL and stat: no need to perform read through
      // 3_read_category/3b_classify_indel
      CLASSIFY_INDEL_NANOPORE ( reads, file(params.ref) )
      ch_versions = ch_versions.mix(CLASSIFY_INDEL_NANOPORE.out.versions)
      STAT_INDEL ( CLASSIFY_INDEL_NANOPORE.out.stat_indel.collect() )
      ch_versions = ch_versions.mix(STAT_INDEL.out.versions)
      // STAT_QC_INDEL ( CLASSIFY_INDEL_NANOPORE.out.reads_no_indel, CLASSIFY_INDEL_NANOPORE.out.mean_qc, CLASSIFY_INDEL_NANOPORE.out.per_site_qc )

      //
      // MODULE: classify_readthrough
      // 3_read_category/3c_classify_readthrough
      // CLASSIFY_READTHROUGH_NANOPORE ( CLASSIFY_INDEL_NANOPORE.out.reads_no_indel )
      // ch_versions = ch_versions.mix(CLASSIFY_READTHROUGH_NANOPORE.out.versions)
      
      // STAT_READTHROUGH ( CLASSIFY_READTHROUGH_NANOPORE.out.stat.collect() )
      // ch_versions = ch_versions.mix(STAT_READTHROUGH.out.versions)
      // FASTQC_READTHROUGH ( reads )
      // ch_versions = ch_versions.mix(FASTQC_READTHROUGH.out.versions)
      // READ_PER_UMI_READTHROUGH ( CLASSIFY_READTHROUGH_NANOPORE.out.reads_through )
      // ch_versions = ch_versions.mix(READ_PER_UMI_READTHROUGH.out.versions)

      // MODULE: classify reads into alleles based on SNP
      // TODO module
      // if params.allele_number == 2:
      //   CLASSIFY_ALLELE ()

    emit:
      reads_no_indel            = CLASSIFY_INDEL_NANOPORE.out.reads_no_indel
      versions      = ch_versions
}

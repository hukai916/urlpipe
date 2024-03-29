include { QC_FLANKING_READS } from '../modules/local/qc_flanking_reads'
include { GET_HIGH_QUALITY_FLANKING_READS } from '../modules/local/get_high_quality_flanking_reads'
include { STAT as STAT_FILTER_BY_INDEL_LENGTH_INDEL_5P_ONLY } from '../modules/local/stat'

include { CLASSIFY_LOCUS              } from '../modules/local/classify_locus'
include { STAT as STAT_LOCUS          } from '../modules/local/stat'

include { CLASSIFY_INDEL_NANOPORE     } from '../modules/local/classify_indel_nanopore'
include { STAT as STAT_INDEL          } from '../modules/local/stat'
include { STAT_QC_INDEL               } from '../modules/local/stat_qc_indel'

include { PARSE_CIGAR as PARSE_CIGAR_INDEL_5P_ONLY } from '../modules/local/parse_cigar'
include { FILTER_BY_INDEL_LENGTH as FILTER_BY_INDEL_LENGTH_INDEL_5P_ONLY  } from '../modules/local/filter_by_indel_length'
include { STAT as STAT_HIGH_QUALITY_FLANKING_READS } from '../modules/local/stat'

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
      // ch_versions = ch_versions.mix (CLASSIFY_LOCUS.out.versions)
      // STAT_LOCUS ( CLASSIFY_LOCUS.out.stat.collect() )
      // ch_versions = ch_versions.mix (STAT_LOCUS.out.versions)

      // MODULE: prefilter to obtain reads with high quality 5' and 3' flanking sequences
      // This module can't be put into preprocessing_nanopore.nf since only one ref is allowed for this module.
      // 3_read_category/3a_high_quality_flanking_read/qc
      QC_FLANKING_READS ( reads, file(params.ref) )
      ch_versions = ch_versions.mix ( QC_FLANKING_READS.out.versions )

      // MODULE: get high quality flanking reads
      // 3_read_category/3a_high_quality_flanking_read
      // QC_FLANKING_READS.out.read_id_mean_qc.view()
      GET_HIGH_QUALITY_FLANKING_READS ( QC_FLANKING_READS.out.reads_input, QC_FLANKING_READS.out.read_id_mean_qc, QC_FLANKING_READS.out.read_id_mean_qc_extend )
      ch_versions = ch_versions.mix ( GET_HIGH_QUALITY_FLANKING_READS.out.versions )
      // 3_read_category/3a_high_quality_flanking_read
      STAT_HIGH_QUALITY_FLANKING_READS ( GET_HIGH_QUALITY_FLANKING_READS.out.stat.collect() )

      //
      // MODULE: classify INDEL and stat: no need to perform read through
      // 3_read_category/3b_classify_indel
      CLASSIFY_INDEL_NANOPORE ( GET_HIGH_QUALITY_FLANKING_READS.out.reads, file(params.ref) )
      ch_versions = ch_versions.mix (CLASSIFY_INDEL_NANOPORE.out.versions)
      STAT_INDEL ( CLASSIFY_INDEL_NANOPORE.out.stat_indel.collect() )
      ch_versions = ch_versions.mix ( STAT_INDEL.out.versions )
      // 3_read_category/3b_classify_indel/qc
      STAT_QC_INDEL ( CLASSIFY_INDEL_NANOPORE.out.reads_no_indel, CLASSIFY_INDEL_NANOPORE.out.mean_qc, CLASSIFY_INDEL_NANOPORE.out.per_site_qc )

      // 
      // MODULE: parse CIGAR from bam and output the ref pos for each CIGAR tuple for each read
      // 3_read_category/3b_classify_indel/indel_5p_only/parse_cigar
      PARSE_CIGAR_INDEL_5P_ONLY ( CLASSIFY_INDEL_NANOPORE.out.reads_indel_5p_only, CLASSIFY_INDEL_NANOPORE.out.bam_minimap2_indel_5p_only )

      // 
      // MODULE: filter indel reads by indel length according to parse_cigar result
      // 3_read_category/3b_classify_indel/indel_5p_only/filter_by_indel_length
      FILTER_BY_INDEL_LENGTH_INDEL_5P_ONLY ( PARSE_CIGAR_INDEL_5P_ONLY.out.reads_input, PARSE_CIGAR_INDEL_5P_ONLY.out.parse_cigar, file(params.ref) )
      // 3_read_category/3b_classify_indel/indel_5p_only/filter_by_indel_length/filter_by_indel_length_indel_5p_only.csv
      STAT_FILTER_BY_INDEL_LENGTH_INDEL_5P_ONLY ( FILTER_BY_INDEL_LENGTH_INDEL_5P_ONLY.out.stat.collect() )

      //
      // MODULE: classify_readthrough
      // 3_read_category/3c_classify_readthrough
      // CLASSIFY_READTHROUGH_NANOPORE ( CLASSIFY_INDEL_NANOPORE.out.reads_no_indel )
      // ch_versions = ch_versions.mix (CLASSIFY_READTHROUGH_NANOPORE.out.versions)
      
      // STAT_READTHROUGH ( CLASSIFY_READTHROUGH_NANOPORE.out.stat.collect() )
      // ch_versions = ch_versions.mix (STAT_READTHROUGH.out.versions)
      // FASTQC_READTHROUGH ( reads )
      // ch_versions = ch_versions.mix (FASTQC_READTHROUGH.out.versions)
      // READ_PER_UMI_READTHROUGH ( CLASSIFY_READTHROUGH_NANOPORE.out.reads_through )
      // ch_versions = ch_versions.mix (READ_PER_UMI_READTHROUGH.out.versions)

      // MODULE: classify reads into alleles based on SNP
      // TODO module
      // if params.allele_number == 2:
      //   CLASSIFY_ALLELE ()

    emit:
      reads_no_indel            = CLASSIFY_INDEL_NANOPORE.out.reads_no_indel
      versions      = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowUrlpipe.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/input_check'
include { USE_READ_R1    } from '../subworkflows/use_read_r1'
include { USE_READ_MERGE } from '../subworkflows/use_read_merge'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { UMI_EXTRACT                 } from '../modules/local/umi_extract'
include { CUTADAPT                    } from '../modules/nf-core/modules/cutadapt/main'
include { FASTQC; FASTQC as FASTQC1   } from '../modules/nf-core/modules/fastqc/main'
include { MAP_LOCUS                   } from '../modules/local/map_locus'
include { CAT_STAT; CAT_STAT as CAT_STAT2; CAT_STAT as CAT_STAT3; CAT_STAT as CAT_STAT5; } from '../modules/local/cat_stat'
include { CAT_STAT_CUTOFF             }   from '../modules/local/cat_stat_cutoff'

include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_2      }   from '../modules/local/cat_stat_cutoff'

include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_INDEL      }   from '../modules/local/cat_stat_cutoff'
include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_INDEL_2}   from '../modules/local/cat_stat_cutoff'

include { UMI_PATTERN; UMI_PATTERN as UMI_PATTERN2 } from '../modules/local/umi_pattern'
include { CLASSIFY_INDEL              } from '../modules/local/classify_indel'
include { CLASSIFY_READTHROUGH        } from '../modules/local/classify_readthrough'

include { FASTQC_SINGLE               } from '../modules/local/fastqc_single'
include { REPEAT_DIST_DISTANCE_R1     } from '../modules/local/repeat_dist_distance_r1'
include { READ_LENGTH_DIST            } from '../modules/local/read_length_dist'

include { REPEAT_DIST_WITHIN_UMI_GROUP as REPEAT_DIST_WITHIN_UMI_GROUP_R1 } from '../modules/local/repeat_dist_within_umi_group'


include { UMI_GROUP_STAT as UMI_GROUP_STAT_R1 } from '../modules/local/umi_group_stat'

include { UMI_GROUP_STAT as UMI_GROUP_STAT_INDEL } from '../modules/local/umi_group_stat'

include { COUNT_SUMMARY as COUNT_SUMMARY_LD                  } from '../modules/local/count_summary'
include { COUNT_SUMMARY as COUNT_SUMMARY_MODE                } from '../modules/local/count_summary'

include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow URLPIPE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input,
        params.allele_number
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Cat Fastq
    //
    CAT_FASTQ (
      INPUT_CHECK.out.reads
      )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    //
    // MODULE: UMI extract
    //
    UMI_EXTRACT (
      CAT_FASTQ.out.reads
      )
    ch_versions = ch_versions.mix(UMI_EXTRACT.out.versions)

    //
    // MODULE: cutadapt
    //
    CUTADAPT (
      UMI_EXTRACT.out.reads
      )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC (
      CUTADAPT.out.reads,
      "0d_fastqc"
      )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    //
    // MODULE: map locus
    //
    MAP_LOCUS (
      CUTADAPT.out.reads
      )
    ch_versions = ch_versions.mix(MAP_LOCUS.out.versions)

    //
    // MODULE: combine MAP_LOCUS.out.stat into one file
    //
    CAT_STAT (
      MAP_LOCUS.out.stat.collect(),
      "1a_map_locus/stat",
      "all_sample",
      "sample_name,locus,percent,misprimed,percent,problem,percent" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT.out.versions)

    //
    // MODULE: UMI pattern: 2a
    //
    UMI_PATTERN (
      CUTADAPT.out.reads,
      "2a_umi_pattern"
      )
    ch_versions = ch_versions.mix(UMI_PATTERN.out.versions)

    //
    // MODULE: classify INDEL
    //
    CLASSIFY_INDEL (
      MAP_LOCUS.out.reads_locus
      )
    ch_versions = ch_versions.mix(CLASSIFY_INDEL.out.versions)

    //
    // MODULE: combine CLASSIFY_INDEL.out.stat into one file
    //
    CAT_STAT2 (
      CLASSIFY_INDEL.out.stat.collect(),
      "3a_classify_indel/stat",
      "all_sample",
      "sample_name,no_indel,no_indel_percent,indel_5p,indel_5p_percent,indel_3p,indel_3p_percent,indel_5p_3p,indel_5p_3p_percent" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT2.out.versions)

    //
    // MODULE: classify_readthrough
    //
    CLASSIFY_READTHROUGH (
      CLASSIFY_INDEL.out.reads_no_indel
      )
    ch_versions = ch_versions.mix(CLASSIFY_READTHROUGH.out.versions)

    //
    // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    //
    CAT_STAT3 (
      CLASSIFY_READTHROUGH.out.stat.collect(),
      "4a_classify_readthrough/stat",
      "all_sample",
      "sample_name,count_readthrough,count_readthrough_percent,count_non_readthrough,p_count_non_readthrough_percent" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT3.out.versions)


    if (params.use_read == "R1") {
      USE_READ_R1( CLASSIFY_READTHROUGH.out.reads_through )
    } else if (params.use_read == "R2") {
      // USE_READ_R2()
      // TODO
      exit 1, '--use_read R2 not implemented yet!'
    } else if (params.use_read == "merge") {
      USE_READ_MERGE( CLASSIFY_READTHROUGH.out.reads_through )
    } else {
      exit 1, '--use_read must be from "R1", "R2", and "merge"!'
    }

 // // For INDEL reads:
    // MODULE: INDEL reads distribution:
    READ_LENGTH_DIST (
      CLASSIFY_INDEL.out.reads_indel_5p_or_3p,
      "4d_indel_read_length_distribution"
    )
    ch_versions = ch_versions.mix(READ_LENGTH_DIST.out.versions)

    // MODULE: INDEL reads UMI stat
    UMI_GROUP_STAT_INDEL (
      READ_LENGTH_DIST.out.count_r1,
      READ_LENGTH_DIST.out.stat_raw, // stat_raw store the raw stat before UMI correction
      "5c_indel_umi_group_stat"
      )
    ch_versions = ch_versions.mix(UMI_GROUP_STAT_INDEL.out.versions)

    // MODULE: UMI correct
    // // if UMI cutoff filter results to 0 counts, the following module hicups, need more debugging
    // REPEAT_DIST_UMI_CORRECT_INDEL (
    //   UMI_GROUP_STAT_INDEL.out.stat,
    //   UMI_GROUP_STAT_INDEL.out.stat_raw,
    //   params.umi_cutoffs,
    //   "5d_indel_read_length_dist_umi_correct"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_UMI_CORRECT_INDEL.out.versions)

    // MODULE: INDEL UMI correct
    READ_UMI_CORRECT (
      UMI_GROUP_STAT_INDEL.out.stat,
      CLASSIFY_INDEL.out.reads_indel_5p_or_3p_pure.collect(),
      params.umi_cutoffs,
      "5e_indel_read_umi_correct"
      )
    ch_versions = ch_versions.mix(READ_UMI_CORRECT.out.versions)

    CAT_STAT_CUTOFF_INDEL ( READ_UMI_CORRECT.out.count_ld.collect(), "ld", "sample_name,read_count",
    "0," + params.umi_cutoffs,
    "all_sample_indel",
    "5e_indel_read_umi_correct/count" )

    CAT_STAT_CUTOFF_INDEL_2 ( READ_UMI_CORRECT.out.count_mode.collect(), "mode", "sample_name,read_count",
    "0," + params.umi_cutoffs,
    "all_sample_indel",
    "5e_indel_read_umi_correct/count" )

    // MODULE: COUNT_SUMMARY: for mode
    COUNT_SUMMARY_MODE (
      USE_READ_R1.out.stat5,
      USE_READ_R1.out.cutoff_stat,
      CAT_STAT_CUTOFF_INDEL_2.out.stat,
      "0," + params.umi_cutoffs,
      "mode",
      "6a_count_summary"
    )
    // MODULE COUNT_SUMMARY: for ld
    COUNT_SUMMARY_LD (
      USE_READ_R1.out.stat5,
      USE_READ_R1.out.cutoff_stat,
      CAT_STAT_CUTOFF_INDEL.out.stat,
      "0," + params.umi_cutoffs,
      "ld",
      "6a_count_summary"
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowUrlpipe.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
    // ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

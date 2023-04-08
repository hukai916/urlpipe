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

include { INPUT_CHECK    } from '../subworkflows/input_check'
include { PREPROCESS_QC  } from '../subworkflows/preprocess_qc'
include { CLASSIFY_READ  } from '../subworkflows/classify_read'

include { REPEAT_STAT_DEFAULT } from '../subworkflows/repeat_stat_default'
include { REPEAT_STAT_MERGE   } from '../subworkflows/repeat_stat_merge'
// include { MODE_NANOPORE  } from '../subworkflows/mode_nanopore'

include { INDEL_STAT     } from '../subworkflows/indel_stat'


include { CAT_STAT; CAT_STAT as CAT_STAT2; } from '../modules/local/cat_stat'
include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_INDEL   }   from '../modules/local/cat_stat_cutoff'
include { CAT_STAT_CUTOFF as CAT_STAT_CUTOFF_INDEL_2 }   from '../modules/local/cat_stat_cutoff'
include { READ_UMI_CORRECT } from '../modules/local/read_umi_correct'
include { READ_LENGTH_DIST } from '../modules/local/read_length_dist'
include { REPEAT_DIST_UMI_CORRECT as REPEAT_DIST_UMI_CORRECT_INDEL } from '../modules/local/repeat_dist_umi_correct'
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
    // 0_pipeline_info/samplesheet.valid.csv
    INPUT_CHECK (
        ch_input,
        params.allele_number
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: preprocess and qc
    // 1_preprocess and 2_qc_and_umi
    PREPROCESS_QC (
      INPUT_CHECK.out.reads,
      ch_versions
      )
    ch_versions = ch_versions.mix(PREPROCESS_QC.out.versions)

    //
    // SUBWORKFLOW: classify reads
    // 3_read_category
    CLASSIFY_READ (
      PREPROCESS_QC.out.reads,
      ch_versions
      )
    ch_versions = ch_versions.mix(CLASSIFY_READ.out.versions)

    // 
    // SUBWORKFLOW: generate repeat statistics
    // 4_repeat_statistics
    if (params.mode == "default") {
      //
      // SUBWORKFLOW: obtain repeat statistics using default mode where individual R1 and R2 reads are used
      // 4_repeat_statistics
      REPEAT_STAT_DEFAULT ( CLASSIFY_READ.out.reads_through, ch_versions )
      ch_versions = REPEAT_STAT_DEFAULT.out.versions
    } else if (params.mode == "merge") {
      // merge mode first merge R1 and R2 reads
      REPEAT_STAT_MERGE ( CLASSIFY_READ.out.reads_through, ch_versions )
      ch_versions = REPEAT_STAT_MERGE.out.versions
    } else if (params.mode == "nanopore") {
      // MODE_NANOPORE ()
      exit 1, '--mode nanopore under development!'
    } else {
      exit 1, '--mode must be from "default", "merge", and "nanopore"!'
    }

    // 
    // SUBWORKFLOW: generate statistics for INDEL reads
    // 5_indel_statistics
    INDEL_STAT ( CLASSIFY_READ.out.reads_indel_5p_or_3p, CLASSIFY_READ.out.reads_indel_5p_or_3p_pure, ch_versions )
    ch_versions = INDEL_STAT.out.versions

    // 1 // MODULE: COUNT_SUMMARY: for mode
    // COUNT_SUMMARY_MODE (
    //   REPEAT_STAT_DEFAULT.out.stat5,
    //   REPEAT_STAT_DEFAULT.out.cutoff_stat,
    //   CAT_STAT_CUTOFF_INDEL_2.out.stat,
    //   "0," + params.umi_cutoffs,
    //   "mode",
    //   "XXX_6a_count_summary"
    // )
    // 2 // MODULE COUNT_SUMMARY: for ld
    // COUNT_SUMMARY_LD (
    //   REPEAT_STAT_DEFAULT.out.stat5,
    //   REPEAT_STAT_DEFAULT.out.cutoff_stat,
    //   CAT_STAT_CUTOFF_INDEL.out.stat,
    //   "0," + params.umi_cutoffs,
    //   "ld",
    //   "XXX_6a_count_summary"
    // )

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

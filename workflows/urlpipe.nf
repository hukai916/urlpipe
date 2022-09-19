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
include { INPUT_CHECK } from '../subworkflows/local/input_check'

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
include { CAT_STAT; CAT_STAT as CAT_STAT2; CAT_STAT as CAT_STAT3; CAT_STAT as CAT_STAT4 } from '../modules/local/cat_stat'
include { UMI_PATTERN                 } from '../modules/local/umi_pattern'
include { CLASSIFY_INDEL              } from '../modules/local/classify_indel'
include { CLASSIFY_READTHROUGH        } from '../modules/local/classify_readthrough'
include { BBMERGE                     } from '../modules/local/bbmerge'
include { FASTQC_SINGLE               } from '../modules/local/fastqc_single'

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
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // INPUT_CHECK.out.reads.view()

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
      "sample_name\tlocus\tpercent\tmisprimed\tpercent\tproblem\tpercent" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT.out.versions)

    //
    // MODULE: UMI pattern
    //
    UMI_PATTERN (
      CUTADAPT.out.reads
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
      "sample_name\tno_indel\tno_indel_percent\tindel_5p\tindel_5p_percent\tindel_3p\tindel_3p_percent\tindel_5p_3p\tindel_5p_3p_percent" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT2.out.versions)

    //
    // MODULE: classify_readthrough
    //
    CLASSIFY_READTHROUGH (
      MAP_LOCUS.out.reads_locus
      )
    ch_versions = ch_versions.mix(CLASSIFY_READTHROUGH.out.versions)

    //
    // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    //
    CAT_STAT3 (
      CLASSIFY_READTHROUGH.out.stat.collect(),
      "4a_classify_readthrough/stat",
      "sample_name\tcount_readthrough\tcount_readthrough_percent\tcount_non_readthrough\tp_count_non_readthrough_percent" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT3.out.versions)

    //
    // MODULE: BBmerge
    //
    BBMERGE (
      CLASSIFY_READTHROUGH.out.reads_through
      )
    ch_versions = ch_versions.mix(BBMERGE.out.versions)

    //
    // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    //
    CAT_STAT4 (
      BBMERGE.out.stat.collect(),
      "4b_bbmerge/stat",
      "sample_name\tcount_merge\tcount_non_merge" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT4.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC1 (
      CLASSIFY_READTHROUGH.out.reads_through,
      "4c_fastqc_r1_r2"
      )
    ch_versions = ch_versions.mix(FASTQC1.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC_SINGLE (
      BBMERGE.out.reads_merged,
      "4c_fastqc_merged"
      )
    ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)

    //
    // MODULE: repeat distribution R1/R2 distance
    //
    REPEAT_DIST_DISTANCE (
      CLASSIFY_READTHROUGH.out.reads_through
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_DISTANCE.out.versions)

    //
    // MODULE: repeat distribution R1 distance
    //
    // REPEAT_DIST_DISTANCE_MERGED (
    //   BBMERGE.out.reads_merged,
    //   "4c_merge_fastqc"
    //   )
    // ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)

    //
    // MODULE: repeat distribution R1 distance
    //
    // REPEAT_DIST_ALIGNMENT (
    //   BBMERGE.out.reads_merged,
    //   "4c_merge_fastqc"
    //   )
    // ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)
    //
    // //
    // // MODULE: repeat distribution R1 distance
    // //
    // REPEAT_DIST_ALIGNMENT_MERGED (
    //   BBMERGE.out.reads_merged,
    //   "4c_merge_fastqc"
    //   )
    // ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)


    // MAP_LOCUS.out.stat.collect()
    //
    // MODULE: umi distribution statistics
    //
    // UMI_STAT (
    //   CUTADAPT.out.reads
    //   )
    // ch_versions = ch_versions.mix(UMI_STAT.out.versions)


    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    //
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
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

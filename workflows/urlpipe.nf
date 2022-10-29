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
include { CAT_STAT; CAT_STAT as CAT_STAT2; CAT_STAT as CAT_STAT3; CAT_STAT as CAT_STAT4; CAT_STAT as CAT_STAT5; CAT_STAT as CAT_STAT6; CAT_STAT as CAT_STAT7; CAT_STAT as CAT_STAT8; CAT_STAT as CAT_STAT9; CAT_STAT as CAT_STAT10; CAT_STAT as CAT_STAT11; CAT_STAT as CAT_STAT5_MERGE; CAT_STAT as CAT_STAT7_MERGE; CAT_STAT as CAT_STAT8_MERGE; CAT_STAT as CAT_STAT9_MERGE; CAT_STAT as CAT_STAT10_MERGE; CAT_STAT as CAT_STAT11_MERGE} from '../modules/local/cat_stat'
include { UMI_PATTERN; UMI_PATTERN as UMI_PATTERN2 } from '../modules/local/umi_pattern'
include { CLASSIFY_INDEL              } from '../modules/local/classify_indel'
include { CLASSIFY_READTHROUGH        } from '../modules/local/classify_readthrough'
include { BBMERGE                     } from '../modules/local/bbmerge'
include { FASTQC_SINGLE               } from '../modules/local/fastqc_single'
include { REPEAT_DIST_DISTANCE        } from '../modules/local/repeat_dist_distance'
include { REPEAT_DIST_DISTANCE_MERGED } from '../modules/local/repeat_dist_distance_merged'
include { REPEAT_DIST_WITHIN_UMI_GROUP as REPEAT_DIST_WITHIN_UMI_GROUP_R1 } from '../modules/local/repeat_dist_within_umi_group'
include { REPEAT_DIST_WITHIN_UMI_GROUP as REPEAT_DIST_WITHIN_UMI_GROUP_MERGE } from '../modules/local/repeat_dist_within_umi_group'

include { UMI_GROUP_STAT as UMI_GROUP_STAT_R1 } from '../modules/local/umi_group_stat'
include { UMI_GROUP_STAT as UMI_GROUP_STAT_MERGE } from '../modules/local/umi_group_stat'
include { REPEAT_DIST_UMI_CORRECT as REPEAT_DIST_UMI_CORRECT_R1 } from '../modules/local/repeat_dist_umi_correct'
include { REPEAT_DIST_UMI_CORRECT as REPEAT_DIST_UMI_CORRECT_MERGE } from '../modules/local/repeat_dist_umi_correct'

include { PLOT_FRAC as PLOT_FRAC_4D_R1    } from '../modules/local/plot_frac'
include { PLOT_FRAC as PLOT_FRAC_4D_MERGE } from '../modules/local/plot_frac'
include { PLOT_FRAC as PLOT_FRAC_5D_R1_FRAC_1    } from '../modules/local/plot_frac'
include { PLOT_FRAC as PLOT_FRAC_5D_R1_FRAC_3    } from '../modules/local/plot_frac'
include { PLOT_FRAC as PLOT_FRAC_5D_R1_FRAC_10   } from '../modules/local/plot_frac'
include { PLOT_FRAC as PLOT_FRAC_5D_R1_FRAC_30   } from '../modules/local/plot_frac'

include { PLOT_FRAC as PLOT_FRAC_5D_MERGE_FRAC_1 } from '../modules/local/plot_frac'
include { PLOT_FRAC as PLOT_FRAC_5D_MERGE_FRAC_3 } from '../modules/local/plot_frac'
include { PLOT_FRAC as PLOT_FRAC_5D_MERGE_FRAC_10} from '../modules/local/plot_frac'
include { PLOT_FRAC as PLOT_FRAC_5D_MERGE_FRAC_30 } from '../modules/local/plot_frac'

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
      "sample_name,count_readthrough,count_readthrough_percent,count_non_readthrough,p_count_non_readthrough_percent" // header to be added
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
      "sample_name,count_merge,count_non_merge" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT4.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC1 (
      CLASSIFY_READTHROUGH.out.reads_through,
      "4c_r1_r2_fastqc"
      )
    ch_versions = ch_versions.mix(FASTQC1.out.versions)

    //
    // MODULE: FastQC
    //
    FASTQC_SINGLE (
      BBMERGE.out.reads_merged,
      "4c_merge_fastqc"
      )
    ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)

    //
    // MODULE: repeat distribution distance with R1/R2 reads
    //
    REPEAT_DIST_DISTANCE (
      CLASSIFY_READTHROUGH.out.reads_through,
      "4d_repeat_distribution_distance"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_DISTANCE.out.versions)

    //
    // MODULE: repeat distribution distance with merged reads
    //
    REPEAT_DIST_DISTANCE_MERGED (
      BBMERGE.out.reads_merged,
      "4d_merge_repeat_distribution_distance"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_DISTANCE_MERGED.out.versions)

    //
    // MODULE: UMI pattern: 5a
    //
    UMI_PATTERN2 (
      CLASSIFY_READTHROUGH.out.reads_through,
      "5a_umi_pattern"
      )
    ch_versions = ch_versions.mix(UMI_PATTERN2.out.versions)


    //
    // MODULE: repeat distribution within umi group
    //
    REPEAT_DIST_WITHIN_UMI_GROUP_R1 (
      REPEAT_DIST_DISTANCE.out.count_r1,
      "5b_r1_repeat_dist_within_umi_group"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_WITHIN_UMI_GROUP_R1.out.versions)

    REPEAT_DIST_WITHIN_UMI_GROUP_MERGE (
      REPEAT_DIST_DISTANCE_MERGED.out.count,
      "5b_merge_repeat_dist_within_umi_group"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_WITHIN_UMI_GROUP_MERGE.out.versions)

    //
    //
    // MODULE: repeat distribution within umi group
    //
    // REPEAT_DIST_WITHIN_UMI_GROUP2 (
    //   REPEAT_DIST_DISTANCE.out.count_r2,
    //   "5b_r2_repeat_dist_within_umi_group"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_WITHIN_UMI_GROUP2.out.versions)


    //
    // MODULE: UMI group stat: UMI read_count mean mode: 5c
    //
    UMI_GROUP_STAT_R1 (
      REPEAT_DIST_DISTANCE.out.count_r1,
      REPEAT_DIST_DISTANCE.out.stat_raw, // stat_raw store the raw stat before UMI correction
      "5c_r1_umi_group_stat"
      )
    ch_versions = ch_versions.mix(UMI_GROUP_STAT_R1.out.versions)

    UMI_GROUP_STAT_MERGE (
      REPEAT_DIST_DISTANCE_MERGED.out.count,
      REPEAT_DIST_DISTANCE_MERGED.out.stat_raw,
      "5c_merge_umi_group_stat"
      )
    ch_versions = ch_versions.mix(UMI_GROUP_STAT_MERGE.out.versions)

    //
    // MODULE: repeat dist UMI corrected: 5d
    //
    // UMI_GROUP_STAT_R1.out.stat.view()
    // UMI_GROUP_STAT_R1.out.stat_raw.view()

    REPEAT_DIST_UMI_CORRECT_R1 (
      UMI_GROUP_STAT_R1.out.stat,
      UMI_GROUP_STAT_R1.out.stat_raw,
      "5d_r1_repeat_dist_umi_correct"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_UMI_CORRECT_R1.out.versions)

    REPEAT_DIST_UMI_CORRECT_MERGE (
      UMI_GROUP_STAT_MERGE.out.stat,
      UMI_GROUP_STAT_MERGE.out.stat_raw,
      "5d_merge_repeat_dist_umi_correct"
      )
    ch_versions = ch_versions.mix(REPEAT_DIST_UMI_CORRECT_MERGE.out.versions)

    //
    // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    //
    // REPEAT_DIST_DISTANCE.out.frac_r1.collect().view()
    CAT_STAT5 (
      REPEAT_DIST_DISTANCE.out.frac_r1.collect(),
      "4d_repeat_distribution_distance/frac_r1",
      "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT5.out.versions)

    //
    // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    //
    CAT_STAT6 (
      REPEAT_DIST_DISTANCE.out.frac_r2.collect(),
      "4d_repeat_distribution_distance/frac_r2",
      "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT6.out.versions)

    CAT_STAT5_MERGE (
      REPEAT_DIST_DISTANCE_MERGED.out.frac.collect(),
      "4d_repeat_distribution_distance/frac_merge",
      "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT5_MERGE.out.versions)

    //
    // MODULE: combine REPEAT_DIST_UMI_CORRECT_R1.out.frac_x into one file
    //
    // REPEAT_DIST_UMI_CORRECT_R1.out.frac_1.collect().view()
    CAT_STAT7 ( REPEAT_DIST_UMI_CORRECT_R1.out.frac_1.collect(), "5d_r1_repeat_dist_umi_correct/frac_1", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT7.out.versions)
    CAT_STAT8 ( REPEAT_DIST_UMI_CORRECT_R1.out.frac_3.collect(), "5d_r1_repeat_dist_umi_correct/frac_3", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT8.out.versions)
    CAT_STAT9 ( REPEAT_DIST_UMI_CORRECT_R1.out.frac_10.collect(), "5d_r1_repeat_dist_umi_correct/frac_10", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT9.out.versions)
    CAT_STAT10 ( REPEAT_DIST_UMI_CORRECT_R1.out.frac_30.collect(), "5d_r1_repeat_dist_umi_correct/frac_30", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT10.out.versions)
    CAT_STAT11 ( REPEAT_DIST_UMI_CORRECT_R1.out.frac_100.collect(), "5d_r1_repeat_dist_umi_correct/frac_100", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT11.out.versions)

    CAT_STAT7_MERGE ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_1.collect(), "5d_merge_repeat_dist_umi_correct/frac_1", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT7.out.versions)
    CAT_STAT8_MERGE ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_3.collect(), "5d_merge_repeat_dist_umi_correct/frac_3", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT8.out.versions)
    CAT_STAT9_MERGE ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_10.collect(), "5d_merge_repeat_dist_umi_correct/frac_10", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT9.out.versions)
    CAT_STAT10_MERGE ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_30.collect(), "5d_merge_repeat_dist_umi_correct/frac_30", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT10.out.versions)
    CAT_STAT11_MERGE ( REPEAT_DIST_UMI_CORRECT_MERGE.out.frac_100.collect(), "5d_merge_repeat_dist_umi_correct/frac_100", "sample_name,blew_count,blew_frac,below_mean,below_std,above_count,above_frac,above_mean,above_std" )
    ch_versions = ch_versions.mix(CAT_STAT11.out.versions)

    // MODULE: plot the mean, std, and frac of all_sample.csv for frac stat.
    PLOT_FRAC_4D_R1 (
      CAT_STAT5.out.stat,
      REPEAT_DIST_UMI_CORRECT_R1.out.stat_raw.collect(),
      "4d_repeat_distribution_distance/plot_r1_frac"
    )

    PLOT_FRAC_4D_MERGE (
      CAT_STAT5_MERGE.out.stat,
      REPEAT_DIST_UMI_CORRECT_MERGE.out.stat_raw.collect(),
      "4d_repeat_distribution_distance/plot_merge_frac"
    )

    PLOT_FRAC_5D_R1_FRAC_1 (
      CAT_STAT7.out.stat,
      REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_1_mode_stat.collect(),
      "5d_r1_repeat_dist_umi_correct/plot_r1_frac_1"
    )
    // PLOT_FRAC_5D_R1_FRAC_3 (
    //   CAT_STAT8.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_3_mode_stat.collect(),
    //   "5d_r1_repeat_dist_umi_correct/plot_r1_frac_3"
    // )
    // PLOT_FRAC_5D_R1_FRAC_10 (
    //   CAT_STAT9.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_10_mode_stat.collect(),
    //   "5d_r1_repeat_dist_umi_correct/plot_r1_frac_10"
    // )
    // PLOT_FRAC_5D_R1_FRAC_30 (
    //   CAT_STAT10.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_30_mode_stat.collect(),
    //   "5d_r1_repeat_dist_umi_correct/plot_r1_frac_30"
    // )
    // PLOT_FRAC_5D_R1_FRAC_100 (
    //   CAT_STAT11.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_R1.out.cutoff_100_mode_stat.collect(),
    //   "5d_r1_repeat_dist_umi_correct/plot_r1_frac_100"
    // )

    // PLOT_FRAC_5D_MERGE_FRAC_1 (
    //   CAT_STAT7_MERGE.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_1_mode_stat.collect(),
    //   "5d_merge_repeat_dist_umi_correct/plot_merge_frac_1"
    // )
    // PLOT_FRAC_5D_MERGE_FRAC_3 (
    //   CAT_STAT8_MERGE.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_3_mode_stat.collect(),
    //   "5d_merge_repeat_dist_umi_correct/plot_merge_frac_3"
    // )
    // PLOT_FRAC_5D_MERGE_FRAC_10 (
    //   CAT_STAT9_MERGE.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_10_mode_stat.collect(),
    //   "5d_merge_repeat_dist_umi_correct/plot_merge_frac_10"
    // )
    // PLOT_FRAC_5D_MERGE_FRAC_30 (
    //   CAT_STAT10_MERGE.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_30_mode_stat.collect(),
    //   "5d_merge_repeat_dist_umi_correct/plot_merge_frac_30"
    // )
    // PLOT_FRAC_5D_MERGE_FRAC_100 (
    //   CAT_STAT11_MERGE.out.stat,
    //   REPEAT_DIST_UMI_CORRECT_MERGE.out.cutoff_100_mode_stat.collect(),
    //   "5d_merge_repeat_dist_umi_correct/plot_merge_frac_100"
    // )


    // MODULE: read length SD plot vs UMI per sample:


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

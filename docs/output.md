# URLpipe: Output
(Some parts adapted from nf-core [TEMPLATE](https://github.com/nf-core/tools/blob/master/nf_core/pipeline-template/docs/output.md).)

## Table of Contents
[Introduction](#introduction)   
[Results](#results)   

## Introduction

This document describes the output produced by URLpipe. By default, all results are saved in the "./results" folder, as specified by the `outdir = "./results"` parameter. Results from different modules are organized into corresponding subdirectories (e.g. "./results/1_preprocess", "./results/2_qc_and_umi", *etc.*)

<!-- todo -->
<!-- A html report summarizing all key results will be generated with MultiQC (via custom plugins) and saved into `./results/multiqc/multiqc_report.html` for quick view. -->

## Results

Nextflow implements a caching mechanism that stores all intermediate and final results in the "./work/" directory. By default, files in the "./results/" are symbolic links to that in the "./work/". To switch from `symlink` to `copy`, use `publish_dir_mode = copy` argument. Below summarizes the main contents of each result folder, using [Example study1](usage.md#example-study1) as an example. Note that since the [sample_dataset1.config](../conf/sample_dataset1.config) specifies `outdir = "./results_dataset1"`, all results will be saved in the ["./results_dataset1"](https://github.com/hukai916/URLpipe_example/tree/main/results_dataset1) directory.

<!-- todo -->
<!-- Also see **module references**: [csv](https://github.com/hukai916/scATACpipe/blob/dev/docs/scATACpipe_module_references.csv) or [xlsx](https://github.com/hukai916/scATACpipe/blob/dev/docs/scATACpipe_module_references.xlsx) for more information. -->

<details markdown="1">
<summary>Subfolder: 0_pipeline_info</summary>

This subfolder contains pipeline execution details and the validated sample sheet file.
* `samplesheet.valid.csv`: Validated samplesheet file.
</details>

<details markdown="1">
<summary>Subfolder: 1_preprocess</summary>

This subfolder contains intermediate files produced during the preprocessing steps.
* `1a_lane_merge/`: Contains merged fastq files for the same libraries.
* `1b_umi_extract/`: Contains fastq files with UMIs extracted and appened to the read names.
* `1c_cutadapt/`: Contains fastq files with adapters removed.
</details>

<details markdown="1">
<summary>Subfolder: 2_qc_and_umi</summary>

This subfolder contains intermediate files generated during the QC and UMI processing steps.
* `2a_fastqc/`: Contains FastQC reports for sequencing data at each preprocessing step.
* `2b_read_per_umi_cutadapt/`: Contains distribution of the number of reads per UMI for fastq files after trimming.
* `2b_read_per_umi_readthrough/`: Contains distribution of the number of reads per UMI for readthrough fastq files.
</details>

<details markdown="1">
<summary>Subfolder: 3_read_category</summary>

This subfolder contains categorized reads and associated statistics.
* `3a_classify_locus/`: Contains reads mapped to the target reference.
  * `classify_locus.csv` ([Exmaple](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/3_read_category/3a_classify_locus/classify_locus.csv)): Statistical summary of on-locus read fractions.
  * `on_target_locus/`: Contains on-locus reads (5'-end of reference appear in R1 and 3'-end of reference appear in R2).
  * `off_target_locus/`: Contains off-locus reads (Neither end of the reference sequence appears in either R1 or R2).
  * `problem_reads/`: Contains problematic reads (One end of the reference sequence appears in either R1 or R2).
* `3b_classify_indel/`: Contains reads categorized by INDEL presence flanking the repeat region.
  * `classify_indel.csv` ([Exmaple](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/3_read_category/3b_classify_indel/classify_indel.csv)): Statistical summarizing of non-INDEL read fractions.
  * `no_indel/`: Contains no-INDEL reads (both sides of repeat flanking sequences appear in either R1 or R2).
  * `indel_5p/`: Contains reads with INDELs occuring in the 5'-end repeat flanking region.
  * `indel_3p/`: Contains reads with INDELs occuring in the 3'-end repeat flanking region.
  * `indel_5p_and_3p/`: Contains reads with INDELs occuring in both the 5'-end and 3'-end repeat flanking regions.
  * `indel_5p_or_3p/`: Contains reads with INDELs occuring in either the 5'-end or 3'-end repeat flanking region.
* `3c_classify_readthrough/`: Contains reads that span the flanking repeat regions at both ends.
  * `classify_readthrough.csv` ([Example](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/3_read_category/3c_classify_readthrough/classify_readthrough.csv)): Statistical summary of readthrough read fractions.
  * `readthrough/`: Contains readthrough reads (both sides of repeat flanking sequences appear in R1).
  * `non_readthrough/`: Contains non-readthrough reads.
  * `stat/`: Contains readthrough read fractions for each individual sample.
</details>

<details markdown="1">
<summary>Subfolder: 4_repeat_statistics</summary>

This subfolder contains repeat length statistics for readthrough reads determined in `3_read_category`.
* `4a_repeat_length_distribution/`: Contains distribution of repeat lengths for each UMI cutoff.
  * `repeat_length_count_default_umi_XXX.csv` ([Example](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/4_repeat_statistics/4a_repeat_length_distribution/repeat_length_count_default_umi_3.csv)): Statistical summary of repeat length distribution when collapsed at each UMI cutoff. 
  *  `repeat_length_count_default_umi_XXX.html` ([Example](https://rawcdn.githack.com/hukai916/URLpipe_example/277d084f233ed41307327bc3b770ecd059a6817b/results_dataset1/4_repeat_statistics/4a_repeat_length_distribution/repeat_length_count_default_umi_3.html)): Bar plot visualizing the repeat length distribution summary.
* `4a_repeat_length_distribution_bwa_length/`: Contains intermediate files used for `length_mode = "reference_align"`.
* `4a_repeat_length_distribution_bwa/`: Contains intermediate BAM and fastq files used for `length_mode = "reference_align"`.
* `4b_repeat_length_distribution_per_umi/`: Contains repeat length distribution per UMI.
  * `csv/` ([Example](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/4_repeat_statistics/4b_repeat_length_distribution_per_umi/csv/repeat_length_distribution_per_umi_3_D103_10uM_R1.csv)): Contains summary tables.
  * `html/` ([Example](https://rawcdn.githack.com/hukai916/URLpipe_example/277d084f233ed41307327bc3b770ecd059a6817b/results_dataset1/4_repeat_statistics/4b_repeat_length_distribution_per_umi/html/repeat_length_distribution_per_umi_3_D103_10uM_R1.html)): Contains bar plots visualizing the above tables.

* `4c_repeat_length_fraction/` ([Example](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/4_repeat_statistics/4c_repeat_length_fraction/repeat_length_fraction_umi_3_D103_10uM_R1.csv)): Contains fraction of repeat lengths falling into different ranges, defined by allele-specific repeat lengths for each sample.
</details>

<details markdown="1">
<summary>Subfolder: 5_indel_statistics</summary>

This subfolder contains repeat length statistics for INDEL reads determined in `3_read_category`.
* `5a_read_count_per_umi_cutoff/` ([Example](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/5_indel_statistics/5a_read_count_per_umi_cutoff/read_count_umi_cutoff_5.csv)): INDEL read count after UMI correction.
</details>

<details markdown="1">
<summary>Subfolder: 6_summary</summary>

This subfolder contains repeat length summary statistics and plots generated by combining results from `4_repeat_statistics`, and `5_indel_statistics`.
* `6a_master_table/`: Statistical tables summarizing repeat lengths per UMI cutoff.
  * `master_table_allele_umi_XXX.csv` ([Example](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/6_summary/6a_master_table/master_table_allele_umi_3.csv)): Fraction of repeat length that falls into different ranges defined with allele-specific repeat length for all samples.
  * `master_table_repeat_bin_umi_XXX.csv` ([Example](https://github.com/hukai916/URLpipe_example/blob/main/results_dataset1/6_summary/6a_master_table/master_table_repeat_bin_umi_3.csv)): Like above but more flexible, stores fraction of repeat length (including INDEL reads) that falls into different bins defined with `repeat_bins = "[(0,50), (51,60), (61,137), (138,154), (155,1000)]"`, for all samples.
* `6b_bin_plot/`: Statistical plots summarizing repeat lengths per UMI cutoff.
  * `master_table_repeat_bin_umi_XXX.count.withoutIndel.html` ([Example](https://rawcdn.githack.com/hukai916/URLpipe_example/277d084f233ed41307327bc3b770ecd059a6817b/results_dataset1/6_summary/6b_bin_plot/master_table_repeat_bin_umi_3.count.withoutIndel.html)): Bin plot for repeat length count for each sample per UMI cutoff excluding INDEL reads.
  * `master_table_repeat_bin_umi_XXX.count.withIndel.html` ([Example](https://rawcdn.githack.com/hukai916/URLpipe_example/277d084f233ed41307327bc3b770ecd059a6817b/results_dataset1/6_summary/6b_bin_plot/master_table_repeat_bin_umi_3.count.withIndel.html)): Bin plot for repeat length count for each sample per UMI cutoff including INDEL reads.
  * `master_table_repeat_bin_umi_XXX.ratio.withoutIndel.html` ([Example](https://rawcdn.githack.com/hukai916/URLpipe_example/277d084f233ed41307327bc3b770ecd059a6817b/results_dataset1/6_summary/6b_bin_plot/master_table_repeat_bin_umi_3.ratio.withoutIndel.html)): Bin plot for repeat length fraction for each sample per UMI cutoff excluding INDEL reads.
  * `master_table_repeat_bin_umi_XXX.ratio.withIndel.html` ([Example](https://rawcdn.githack.com/hukai916/URLpipe_example/277d084f233ed41307327bc3b770ecd059a6817b/results_dataset1/6_summary/6b_bin_plot/master_table_repeat_bin_umi_3.ratio.withIndel.html)): Bin plot for repeat length fraction for each sample per UMI cutoff including INDEL reads.
</details>
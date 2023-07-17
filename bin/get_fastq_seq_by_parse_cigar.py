#!/usr/bin/env python

"""
Obtain fastq read given: input fastq, parse_cigar_csv, indel_length_cutff, ref_range, out_pass_reads.fastq.gz, out_not_pass_reads.fastq.gz

Usage:
    python get_fastq_seq_by_parse_cigar.py sample.fastq.gz parse_cigar.csv, indel_length_cutof, ref_range, out_pass.fastq.gz, out_not_pass.fastq.gz 
"""

import sys
from utils import _open, _open_out
from Bio import SeqIO

reads = sys.argv[1]
parse_cigar = sys.argv[2]
indel_length_cutoff = int(sys.argv[3])
ref_start, ref_end = [int(x) for x in str(sys.argv[4]).strip().split(":")]
out_pass_reads = sys.argv[5]
out_not_pass_reads = sys.argv[6]

# Step1: figure out which reads to keep using indel_length_cutoff and ref_range
filter_fastq = {}
with open(parse_cigar, "r") as f:
    for line in f:
        tem = line.strip().split(",")
        read_id, cigar_1, cigar_2, ref_pos = tem[0], tem[1], int(tem[2]), int(tem[3])
        if cigar_1 in ["D", "I"] and cigar_2 >= indel_length_cutoff and ref_pos >= ref_start and ref_pos < ref_end:
            filter_fastq.append(read_id)

# Step2: output to corresponding fastq files
with _open(reads) as f, _open_out(out_pass_reads) as f_out_pass, _open_out(out_not_pass_reads) as f_out_not_pass:
    for record in SeqIO.parse(f, "fastq"):
        if record.id in filter_fastq:
            SeqIO.write(record, f_out_pass)
        else:
            SeqIO.write(record, f_out_not_pass)
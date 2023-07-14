#!/usr/bin/env python

"""
Obtain fastq read given: input fastq, parse_cigar_csv, indel_length_cutff, ref_range, out_pass_reads.fastq.gz, out_not_pass_reads.fastq.gz

Usage:
    python get_fastq_seq_by_parse_cigar.py sample.fastq.gz parse_cigar.csv, indel_length_cutof, ref_range, out_pass.fastq.gz, out_not_pass.fastq.gz 
"""

import sys
from utils import _open
from Bio import SeqIO

reads = sys.argv[1]
parse_cigar = sys.argv[2]
indel_length_cutoff = int(sys.argv[3])
ref_start, ref_end = str(sys.argv[4]).strip().split(":")
out_pass_reads = sys.argv[5]
out_not_pass_reads = sys.argv[6]

# Step1: figure out which reads to keep using indel_length_cutoff and ref_range
filter_fastq = {}
with open(parse_cigar, "r") as f:
    for line in f:
        tem = line.strip().split(",")



range_dict = {}

with open(start_pos, "r") as f:
    for line in f:
        tem = line.strip().split(",")
        if tem[0] not in range_dict:
            range_dict[tem[0]] = [tem[1]]
with open(end_pos, "r") as f:
    for line in f:
        tem = line.strip().split(",")
        if tem[0] in range_dict:
            range_dict[tem[0]].append(tem[1])

with _open(reads) as f, open(outfile, "w") as f_out:
    for record in SeqIO.parse(f, "fastq"):
        if record.id in range_dict:
            if len(range_dict[record.id]) != 2: # only start exist, but not end
                f_out.write(record.id + ",NA\n")
                # print(record.id, range_dict[record.id])
            elif range_dict[record.id][0] == "NA" or range_dict[record.id][1] == "NA":
                f_out.write(record.id + ",NA\n")
            else:
                # qc_score = [str(x) for x in record.letter_annotations["phred_quality"][int(range_dict[record.id][0]):int(range_dict[record.id][1])]]
                seq = [str(x) for x in record.seq[int(range_dict[record.id][0]):int(range_dict[record.id][1])]]
                
                f_out.write(record.id + "," + "".join(seq) + "\n")


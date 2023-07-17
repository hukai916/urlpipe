#!/usr/bin/env python

"""
Return reads with high quality flanking regions.

Usage:

python get_high_quality_flanking_reads.py sample_reads.fastq.gz read_id_mean_qc mean_qc_cutoff pass_name.fastq.gz not_pass_name.fastq.gz mode read_id_mean_qc_extend

the "mode" can be either "strict" or "default"; if "strict", will also check read_id_mean_qc_extend file, which extends 10nt downstream for left_flanking reads and 10nt upstream for right_flanking reads, and see if mean quality of them also above cutoff.


"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import gzip
from mimetypes import guess_type
from functools import partial
from utils import _open, _open_out

fastq = sys.argv[1]
read_id_mean_qc = sys.argv[2]
mean_qc_cutoff = int(sys.argv[3])
out_pass_name = sys.argv[4]
out_not_pass_name = sys.argv[5]
mode = sys.argv[6]
read_id_mean_qc_extend = sys.argv[7]

# step1: read in read_id_mean_qc
read_id_mean_qc_dict = {}
with open(read_id_mean_qc, "r") as f:
    for line in f:
        tem = line.strip().split(",")
        read_id_mean_qc_dict[tem[0]] = [float(tem[1]), float(tem[2])]

# step2: read in read_id_mean_qc_extend 
read_id_mean_qc_extend_dict = {}
with open(read_id_mean_qc_extend, "r") as f:
    for line in f:
        tem = line.strip().split(",")
        read_id_mean_qc_extend_dict[tem[0]] = [float(tem[1]), float(tem[2])]

# step2: filter reads based on mean_qc
with _open(fastq) as f_in, _open_out(out_pass_name) as f_out_pass, _open_out(out_not_pass_name) as f_out_not_pass:
    for record in SeqIO.parse(f_in, "fastq"):
        if record.id in read_id_mean_qc_dict:
            if mode == "default":
                if read_id_mean_qc_dict[record.id][0] >= mean_qc_cutoff and read_id_mean_qc_dict[record.id][1] >= mean_qc_cutoff:
                    SeqIO.write(record, f_out_pass, "fastq")
                else:
                    SeqIO.write(record, f_out_not_pass, "fastq")
            elif mode == "strict":
                if read_id_mean_qc_dict[record.id][0] >= mean_qc_cutoff and read_id_mean_qc_dict[record.id][1] >= mean_qc_cutoff and read_id_mean_qc_extend_dict[record.id][0] >= mean_qc_cutoff and read_id_mean_qc_extend_dict[record.id][1] >= mean_qc_cutoff:
                    SeqIO.write(record, f_out_pass, "fastq")
                else:
                    SeqIO.write(record, f_out_not_pass, "fastq")
        else:
            SeqIO.write(record, f_out_not_pass, "fastq")
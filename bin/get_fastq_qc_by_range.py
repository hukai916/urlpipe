#!/usr/bin/env python

"""
Obtain mean read quality by range.

Usage:
    python get_fastq_mean_qc_by_range.py sample.fastq.gz start_pos.txt end_pos.txt out_qc.txt
"""

import sys
from utils import _open
from Bio import SeqIO

reads = sys.argv[1]
start_pos = sys.argv[2]
end_pos = sys.argv[3]
outfile = sys.argv[4]

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
                qc_score = [str(x) for x in record.letter_annotations["phred_quality"][int(range_dict[record.id][0]):int(range_dict[record.id][1])]]
                f_out.write(record.id + "," + ",".join(qc_score) + "\n")


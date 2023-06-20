#!/usr/bin/env python

"""
Return base in given range.

Usage:

python classify_indel_nanopore_location_info.py \\
        sample_reads.fastq.gz sample_name \\
        repeat_flanking_left_in_range.txt \\
        repeat_flanking_right_in_range.txt \\ 
        no_indel/${prefix}.fastq.gz \\
        indel_5p_only/${prefix}.fastq.gz \\
        indel_3p_only/${prefix}.fastq.gz \\
        indel_5p_and_3p/${prefix}.fastq.gz \\
        undetermined/${prefix}.fastq.gz \\
        stat/${prefix}.csv
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import gzip
from mimetypes import guess_type
from functools import partial

fastq = sys.argv[1]
sample_name = sys.argv[2]
repeat_flanking_left = sys.argv[3]
repeat_flanking_right = sys.argv[4]
out_no_indel = sys.argv[5]
out_indel_5p_only = sys.argv[6]
out_indel_3p_only = sys.argv[7]
out_indel_5p_and_3p = sys.argv[8]
out_undetermined = sys.argv[9]
out_stat = sys.argv[10]

encoding = guess_type(fastq)[1]
encoding_out = guess_type(out_no_indel)[1]

_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
_open_out = partial(gzip.open, mode='wt') if encoding_out == 'gzip' else open

# obtain matched reads id
match_dict = {}
for line in open(repeat_flanking_left):
    tem = line.split()
    if not tem[1] == "NA":
        match_dict[tem[0]] = [int(tem[1]), "NA"]
for line in open(repeat_flanking_right):
    tem = line.split()
    if not tem[1] == "NA":
        if tem[0] in match_dict:
            match_dict[tem[0]][1] = int(tem[1])
        else:
            match_dict[tem[0]] = ["NA", int(tem[1])]

# output reads
total_count = 0
no_indel_count = 0
indel_5p_only_count = 0
indel_3p_only_count = 0
indel_5p_and_3p_count = 0
undetermined_count = 0
with _open(fastq) as f_in, _open_out(out_no_indel) as f_out_no_indel, _open_out(out_indel_5p_only) as f_out_indel_5p_only, _open_out(out_indel_3p_only) as f_out_indel_3p_only, _open_out(out_indel_5p_and_3p) as f_out_indel_5p_and_3p, _open_out(out_undetermined) as f_out_undetermined:
    # not appear: NA and NA: undetermined
    # NA int
    # int NA 
    # int int 
    #     if int1 > int2
    #         undetermined 
    #     else
    #         no_indel 
    
    for record in SeqIO.parse(f_in, "fastq"):
        total_count += 1
        if not record.id in match_dict:
            indel_5p_and_3p_count += 1
            SeqIO.write(record, f_out_indel_5p_and_3p, "fastq")
        elif match_dict[record.id][0] == "NA":
            indel_5p_only_count += 1
            SeqIO.write(record, f_out_indel_5p_only, "fastq")
        elif match_dict[record.id][1] == "NA":
            indel_3p_only_count += 1
            SeqIO.write(record, f_out_indel_3p_only, "fastq")
        else:
            if match_dict[record.id][0] > match_dict[record.id][1]:
                undetermined_count += 1
                SeqIO.write(record, f_out_undetermined, "fastq")
            else:
                no_indel_count += 1
                SeqIO.write(record, f_out_no_indel, "fastq")

with open(out_stat, "w") as f:
    res = [sample_name]
    res.append(str(no_indel_count))
    res.append(str(indel_5p_only_count))
    res.append(str(indel_3p_only_count))
    res.append(str(indel_5p_and_3p_count))
    res.append(str(undetermined_count))
    res.append(str(no_indel_count / total_count))
    res.append(str(indel_5p_only_count / total_count))
    res.append(str(indel_3p_only_count / total_count))
    res.append(str(indel_5p_and_3p_count / total_count))
    res.append(str(undetermined_count / total_count))
    
    f.write(",".join(res) + "\n")
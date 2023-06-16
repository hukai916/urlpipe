#!/usr/bin/env python

"""
Return base in given range.

Usage:

python get_full_length_reads.py sample_reads.fastq.gz ref_start_in_range.txt ref_end_in_range.txt read_start_range read_end_range


"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import gzip
from mimetypes import guess_type
from functools import partial

fastq = sys.argv[1]
ref_start_in_range = sys.argv[2]
ref_end_in_range = sys.argv[3]
read_start_range = sys.argv[4]
read_end_range = sys.argv[5]
out_in_range_reads = sys.argv[6]
out_not_in_range_reads = sys.argv[7]

encoding = guess_type(fastq)[1]
encoding_out = guess_type(out_in_range_reads)[1]

_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
_open_out = partial(gzip.open, mode='wt') if encoding_out == 'gzip' else open


# obtain: read_start_s, read_start_e, read_end_s, read_end_e
_tem_s = read_start_range.split(":")
if _tem_s[0] == "":
    read_start_s = 0
else:
    read_start_s = int(_tem_s[0])
read_start_e = int(_tem_s[1])

_tem_s = read_end_range.split(":")
if _tem_s[1] == "":
    read_end_e = 0
else:
    read_end_e = int(_tem_s[1])
read_end_s = int(_tem_s[0])

# print(read_start_s, read_start_e, read_end_s, read_end_e)

# obtain in_range reads id
ref_start_range_dict = {}
for line in open(ref_start_in_range):
    tem = line.split()
    if not tem[1] == "NA":
        if read_start_s <= int(tem[1]) <= read_start_e:
            ref_start_range_dict[tem[0]] = 1

ref_end_range_dict = {}
for line in open(ref_end_in_range):
    tem = line.split()
    if not tem[1] == "NA":
        if read_end_s <= int(tem[1]) - int(tem[2]) <= read_end_e:
            ref_end_range_dict[tem[0]] = 1
            
read_in_range = {}
for key in ref_start_range_dict:
    if key in ref_end_range_dict:
        read_in_range[key] = 1

# output in_range and not_in_range reads
with _open(fastq) as f_in, _open_out(out_in_range_reads) as f_out_in_range, _open_out(out_not_in_range_reads) as f_out_not_in_range:
    for record in SeqIO.parse(f_in, "fastq"):
        if record.id in read_in_range:
            SeqIO.write(record, f_out_in_range, "fastq")
        else:
            SeqIO.write(record, f_out_not_in_range, "fastq")
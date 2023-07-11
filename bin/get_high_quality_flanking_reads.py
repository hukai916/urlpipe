#!/usr/bin/env python

"""
Return reads with high quality flanking regions.

Usage:

python get_high_quality_flanking_reads.py sample_reads.fastq.gz read_id_mean_qc mean_qc_cutoff pass_name.fastq.gz not_pass_name.fastq.gz


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

_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
_open_out = partial(gzip.open, mode='wt') if encoding_out == 'gzip' else open

# step1: read in read_id_mean_qc
read_id_mean_qc_dict = {}
with open(read_id_mean_qc, "r") as f:
    for line in f:
        tem = line.strip().split(",")
        read_id_mean_qc_dict[tem[0]] = int(tem[1])

# step2: filter reads based on mean_qc
with _open(fastq) as f_in, _open_out(out_pass_name) as f_out_pass, _open_out(not_pass_name) as f_out_not_pass:
    for record in SeqIO.parse(f_in, "fastq"):
        if record.id in read_id_mean_qc_dict:
            if read_id_mean_qc_dict[record.id] >= mean_qc_cutoff:
                SeqIO.write(record, f_out_pass, "fastq")
            else:
                SeqIO.write(record, f_out_not_pass, "fastq")
        else:
            SeqIO.write(record, f_out_not_pass, "fastq")
#!/usr/bin/env python

"""
Split reads into snp1, snp2, and undetermined using snp information.

Usage:

python split_allele.py sample_reads.fastq.gz snp1.csv snp2.csv split_allele/snp1_$reads split_allele/snp2_$reads split_allele/undetermined_$reads

"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import gzip
from mimetypes import guess_type
from functools import partial

fastq = sys.argv[1]
snp1  = sys.argv[2]
snp2  = sys.argv[3]
out_snp1 = sys.argv[4]
out_snp2 = sys.argv[5]
out_undetermined = sys.argv[6]

encoding = guess_type(fastq)[1]
encoding_out = guess_type(out_snp1)[1]

_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
_open_out = partial(gzip.open, mode='wt') if encoding_out == 'gzip' else open

# obtain snp1/snp2 reads id
snp1_dict = {}
for line in open(snp1):
    tem = line.split()
    if not tem[1] == "NA":
        snp1_dict[tem[0]] = 1
snp2_dict = {}
for line in open(snp2):
    tem = line.split()
    if not tem[1] == "NA":
        snp2_dict[tem[0]] = 1

# output snp1, snp2, and undetermined reads
# undetermined: if both snp1 and snp2, or if neither snp1 nor snp2
with _open(fastq) as f_in, _open_out(out_snp1) as f_out_snp1, _open_out(out_snp2) as f_out_snp2, _open_out(out_undetermined) as f_out_undetermined:
    for record in SeqIO.parse(f_in, "fastq"):
        if record.id in snp1_dict and record.id in snp2_dict:
            SeqIO.write(record, f_out_undetermined, "fastq")
        elif record.id in snp1_dict:
            SeqIO.write(record, f_out_snp1, "fastq")
        elif record.id in snp2.dict:
            SeqIO.write(record, f_out_snp2, "fastq")
        else:
            SeqIO.write(record, f_out_undetermined, "fastq")
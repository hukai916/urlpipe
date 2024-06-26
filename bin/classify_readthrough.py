#!/usr/bin/env python

"""
To classify reads into read through and non-read-through.
As long as R1 or R2 contain both repeat flanking reads, the read pair will be treated as read-through read pair.

Usage:
    classify_readthrough.py readthrough non_readthrough stat sample_name $args

"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os
import regex
import gzip
from mimetypes import guess_type
from functools import partial
import regex
import csv

r1 = sys.argv[1]
r2 = sys.argv[2]
readthrough_dir = sys.argv[3]
non_readthrough_dir = sys.argv[4]
stat_dir = sys.argv[5]
sample_name = sys.argv[6]
r1_flanking = sys.argv[7]
r2_flanking = sys.argv[8]
mismatch    = sys.argv[9]

r1_readthrough = {}
r2_readthrough = {}
r_match  = {}

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

r1_dna = Seq(r1_flanking)
r2_dna = Seq(r2_flanking)
r1_flanking_rc = str(r1_dna.reverse_complement())
r2_flanking_rc = str(r2_dna.reverse_complement())

with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r1_flanking + ")" + "{e<=" + str(mismatch) + "}", str(record.seq)) and regex.search("(" + r2_flanking + ")" + "{e<=" + str(mismatch) + "}", str(record.seq)):
            r1_readthrough[record.name] = 1
        else:
            r1_readthrough[record.name] = 0
with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if regex.search("(" + r2_flanking_rc + ")" + "{e<=" + str(mismatch) + "}", str(record.seq)) and regex.search("(" + r1_flanking_rc + ")" + "{e<=" + str(mismatch) + "}", str(record.seq)):
            r2_readthrough[record.name] = 1
        else:
            r2_readthrough[record.name] = 0

for k in r1_readthrough:
    assert k in r2_readthrough, "not equal records!"
    r_match[k] = r1_readthrough[k] + r2_readthrough[k]

# output:
out_readthrough_r1 = os.path.join(readthrough_dir, sample_name + "_1.fastq")
out_readthrough_r2 = os.path.join(readthrough_dir, sample_name + "_2.fastq")
out_non_readthrough_r1 = os.path.join(non_readthrough_dir, sample_name + "_1.fastq")
out_non_readthrough_r2 = os.path.join(non_readthrough_dir, sample_name + "_2.fastq")
stat_outfile = os.path.join(stat_dir, sample_name + ".csv")

for dir in [out_readthrough_r1, out_readthrough_r2, out_non_readthrough_r1, out_non_readthrough_r2, stat_outfile]:
    os.makedirs(os.path.dirname(dir), exist_ok = True)

r_readthrough_set = set()
r_non_readthrough_set = set()
with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] >= 1:
            r_readthrough_set.add(record.name)
        else:
            r_non_readthrough_set.add(record.name)
with _open(r2) as f:
    for record in SeqIO.parse(f, 'fastq'):
        if r_match[record.name] >= 1:
            r_readthrough_set.add(record.name)
        else:
            r_non_readthrough_set.add(record.name)

with open(out_readthrough_r1, "w") as output_r1_readthrough, \
     open(out_non_readthrough_r1, "w") as output_r1_non_readthrough, \
     open(out_readthrough_r2, "w") as output_r2_readthrough, \
     open(out_non_readthrough_r2, "w") as output_r2_non_readthrough:

    with _open(r1) as f1, _open(r2) as f2:
        for record_r1, record_r2 in zip(SeqIO.parse(f1, 'fastq'), SeqIO.parse(f2, 'fastq')):
            if record_r1.name in r_readthrough_set:
                output_r1_readthrough.write(record_r1.format("fastq"))
                output_r2_readthrough.write(record_r2.format("fastq"))
            else:
                output_r1_non_readthrough.write(record_r1.format("fastq"))
                output_r2_non_readthrough.write(record_r2.format("fastq"))

# print some stats:
with open(stat_outfile, "w", newline = '') as f:
    writer = csv.writer(f)
    count_through = len(r_readthrough_set)
    count_non_through = len(r_non_readthrough_set)
    p_count_through = count_through/(sum([count_through, count_non_through]) + 0.1)
    p_count_non_through = count_non_through/(sum([count_through, count_non_through]) + 0.1)

    writer.writerow([sample_name, count_through, count_non_through, p_count_through, p_count_non_through])

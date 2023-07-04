#!/usr/bin/env python

"""
Get mean quality per read.

Usage:
python get_mean_qc.py sample.fastq.gz sample_X output.txt

"""

from Bio import SeqIO
import sys
import gzip
from mimetypes import guess_type
from functools import partial

fastq = sys.argv[1]
sample_name = sys.argv[2]
outfile = sys.argv[3]

encoding = guess_type(fastq)[1]
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

mean_qualities = []

with _open(fastq) as f:
    for record in SeqIO.parse(f, "fastq"):
        mean_quality = sum(record.letter_annotations["phred_quality"]) / len(record)
        mean_qualities.append(mean_quality)

with open(outfile, "w") as f:
    res = ",".join([str(x) for x in mean_qualities])
    f.write(sample_name + "," + res + "\n")
    
#!/usr/bin/env python

"""
output: UMI,read count
"""

import sys
from Bio import SeqIO
import sys
import gzip
from mimetypes import guess_type
from functools import partial

r1 = sys.argv[1]
outfile = sys.argv[2]

encoding = guess_type(r1)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

d_umi = {}
with _open(r1) as f:
    for record in SeqIO.parse(f, 'fastq'):
        umi = record.name.split("_")[1]
        if not umi in d_umi:
            d_umi[umi] = 1
        else:
            d_umi[umi] = d_umi[umi] + 1

with open(outfile, "w") as f:
    for umi in d_umi:
        f.write(umi + "," + str(d_umi[umi]) + "\n")
